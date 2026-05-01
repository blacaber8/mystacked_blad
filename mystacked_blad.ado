
capture mata: mata drop mahal_dist()

mata:
real colvector mahal_dist(string scalar dlist, string scalar Sinvname)
{
    real matrix X, Sinv
    real colvector md
    real scalar i, n

    Sinv = st_matrix(Sinvname)
    X    = st_data(., tokens(dlist))
    n    = rows(X)
    md   = J(n,1,.)

    for (i=1; i<=n; i++) {
        md[i] = sqrt(X[i,.] * Sinv * X[i,.]')
    }

    return(md)
}
end


cap program drop mystacked_blad 
program mystacked_blad, eclass 

syntax    [if] [in] , ///
		  rperiod(integer) event(varname numeric)  xvar(varlist) time(varname numeric) id(varname)  ///
		  [ force    path_actual  by(varlist)  seed(numlist integer max=1) noyet(numlist integer max=1) allstack onlynoyet ] 
		 

quietly {
preserve 
marksample touse
keep if `touse'


* Se o usuário não informar, usa pasta temporária do Stata
local open=0 
if "`path_actual'" == "" {
    local tempdir `"`c(tmpdir)'"'
	local open = 1 
}

if "`path_actual'" != "" {
local tempdir `"`c(pwd)'"'
local open = 0 
}

if "`by'"== "" {
gen _FS=1
local by "_FS"    
}

display "Temporary files will be stored in: `tempdir'"
tempvar tag_t tag_c nt nc prod gtag
egen `tag_t' = tag(`by' `id' `event') if `event' > 0
egen `tag_c' = tag(`by' `id' `event') if `event' == 0

bysort `by': egen `nt' = total(`tag_t')
bysort `by': egen `nc' = total(`tag_c')

gen double `prod' = `nt' * `nc'
egen `gtag' = tag(`by')

quietly summarize `prod' if `gtag', meanonly
local c = r(sum)

dis as red  "Expected treated-control combinations: `c'"
if `c' > 50000000 &  "`force'" == "" {
    di as err "Execution stopped: the number of treated-control combinations (`c') exceeds the allowed limit of 50,000,000."
    di as err "Number of treated units: `a'"
    di as err "Number of control units: `b'"
    di as err "Use a more restrictive sample or set the force option to proceed."
	restore 
    exit 198
}

capture isid `id' `time' `event'
if _rc>0 {
    di as err "The data are not uniquely identified by id and time."
    di as err "There are multiple observations within at least one id-time cell."
	restore 
    exit 459
}

**** Execution:  

* Baseline sample 
saveold  "`tempdir'\\baseline_sample.dta", replace 

* Matrix for mahalanobis 
use "`tempdir'/baseline_sample.dta", clear
keep if `event' == 0
quietly corr `xvar', covariance
matrix S = r(C)

capture matrix Sinv = invsym(S)
if _rc!=0  {
    di as err "Mahalanobis distance could not be computed because the covariance matrix of xvar() is singular."
    di as err "Try removing collinear or near-constant covariates."
	restore 
    exit 498
}

* Control sample 
use "`tempdir'\\baseline_sample.dta", clear 

if "`noyet'"!="" {
keep `id' `by' `event'
bysort `id' `by': keep if _n==1 
ren `id' `id'0 
ren `event' `event'0 
gen ID=1 	
}

if "`noyet'"=="" {
keep if `event'==0 
keep `id' `by'
bysort `id' `by': keep if _n==1 
ren `id' `id'0 
gen ID=1 
}

saveold "`tempdir'\\control_group.dta", replace 

* Treatment sample 
use "`tempdir'\\baseline_sample.dta", clear 
keep if `event'>0 
keep `id' `event' `by'
ren `id' `id'1
bysort `id'1 `event' `by': keep if _n==1 

gen ID=1 

joinby ID `by' using "`tempdir'\\control_group.dta" 
erase  "`tempdir'\\control_group.dta" 

if "`noyet'"!="" {
drop if `id'0 ==`id'1
keep if `event'0-`event'>=`noyet' | `event'0==0  
}

* Only full stacked
if "`allstack'"!="" {

if "`onlynoyet'"!="" {
gen noyet=`event'0>0
}

if "`onlynoyet'"=="" {
gen noyet=1
}
	
keep `id'* `event' noyet
gen match_pair=_n 
gen stacked=`id'1

reshape long  `id' , i(match_pair `event' stacked noyet)  j(treated)
drop match_pair 
ren stacked match_pair 

if "`onlynoyet'"!="" {
keep if noyet==1
drop noyet 
}

if "`onlynoyet'"=="" {
drop noyet 
}

bysort match_pair `id' `event': keep if _n==1
bysort `id' `event': gen weight=_N 
bysort `id' `event': keep if _n==1 

label var weight "Weight for estimation"
label var match_pair "New pair id code"
label var `id' "original ID" 
label var `event' "original event" 
label var treated "New treated/control status"

saveold "`tempdir'\\final_stacked_matched.dta", replace 
erase "`tempdir'\\baseline_sample.dta" 
restore  

di as red "Matching completed successfully."
	
use "`tempdir'\\final_stacked_matched.dta", clear 
erase "`tempdir'\\final_stacked_matched.dta" 
exit
}	

drop `event'0
* Matched on covariates 
gen `time'=`event'-`rperiod'
drop ID 

* Merging covariates 
foreach i in 1 0 {
ren `id'`i' `id'    

if "`i'"=="0" {
	ren `event' _X
	gen `event'=0
}

merge m:1 `id' `time' `event' using "`tempdir'\\baseline_sample.dta", keep(3) keepusing(`xvar') nogen 


if "`i'"=="0" {
	drop `event'
	ren _X `event' 
}

foreach j in  `xvar' {
ren `j' `j'`i'    
}

ren `id' `id'`i'      
}
 
* Mahalanobis distance 
local dlist
local j = 0
foreach v in `xvar' {
    local ++j
    gen d`j' = `v'0 - `v'1
    local dlist `dlist' d`j'
}

gen double maha = .
mata: st_store(., "maha", mahal_dist(st_local("dlist"), "Sinv"))

if "`seed'"!="" {
	set seed `seed'
}

if "`seed'"=="" {
	set seed 123
}	
	gen random=runiform()
	local random "random"


count if maha==. 
if r(N)>0 {
	dis as red "Maha score has missing values"
}

* nearest neighbor 
sort `id'1 `event' maha  `random'
by `id'1 `event': keep if _n==1 

keep `id'* `event' 
gen match_pair=_n 

reshape long  `id' , i(match_pair `event')  j(treated)

label var match_pair "New pair id code"
label var `id' "original ID" 
label var `event' "original event" 
label var treated "New treated/control status"

saveold "`tempdir'\\final_stacked_matched.dta", replace 
erase "`tempdir'\\baseline_sample.dta" 
restore  

di as red "Matching completed successfully."
	
use "`tempdir'\\final_stacked_matched.dta", clear 
erase "`tempdir'\\final_stacked_matched.dta" 

}
describe *

end 
