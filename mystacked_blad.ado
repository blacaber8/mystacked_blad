
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
		  [ force    path_actual  by(varlist)  seed(numlist integer max=1) ] 
		 

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
egen `tag_t' = tag(`by' `id') if `event' > 0
egen `tag_c' = tag(`by' `id') if `event' == 0

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

capture isid `id' `time'
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

use "`tempdir'\\baseline_sample.dta", clear 
* Control sample 
keep if `event'==0 

keep `id' `by'
bysort `id': keep if _n==1 
ren `id' `id'0 
gen ID=1 
saveold "`tempdir'\\control_group.dta", replace 

* Treatment sample 
use "`tempdir'\\baseline_sample.dta", clear 
keep if `event'>0 
keep `id' `event' `by'
ren `id' `id'1
bysort `id'1 `event': keep if _n==1 

gen ID=1 

joinby ID `by' using "`tempdir'\\control_group.dta" 
erase  "`tempdir'\\control_group.dta" 

gen `time'=`event'-`rperiod'
drop ID 

* Merging covariates 
foreach i in 1 0 {
ren `id'`i' `id'    

merge m:1 `id' `time'  using "`tempdir'\\baseline_sample.dta", keep(3) keepusing(`xvar') nogen 

foreach j in  `xvar' {
ren `j' `j'`i'    
}

ren `id' `id'`i'      
}

/*
* computing Euclidian distance 
gen maha=0 
foreach i in `xvar' {
 bysort `id'1: egen SD=sd(`i'0)
 replace maha=maha+ ((`i'0-`i'1)/SD)^2 
 drop SD 
}
replace maha =sqrt(maha)
*/ 
 
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
	gen random=runiform()
	local random "random"
}

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

/*
clear
set obs 20

* 5 unidades, 4 periodos cada = 20 observacoes
gen id = ceil(_n/4)
bysort id: gen time = _n

* Evento: ids 1 e 2 tratados; ids 3, 4 e 5 controles
gen event = .
replace event = 3 if id == 1
replace event = 2 if id == 2
replace event = 0 if inlist(id, 3, 4, 5)

* Covariaveis observadas ao longo do tempo
gen x1 = .
gen x2 = .

replace x1 = 10 if id == 1 & time == 1
replace x1 = 12 if id == 1 & time == 2
replace x1 = 14 if id == 1 & time == 3
replace x1 = 15 if id == 1 & time == 4

replace x1 = 20 if id == 2 & time == 1
replace x1 = 21 if id == 2 & time == 2
replace x1 = 23 if id == 2 & time == 3
replace x1 = 24 if id == 2 & time == 4

replace x1 = 11 if id == 3 & time == 1
replace x1 = 13 if id == 3 & time == 2
replace x1 = 15 if id == 3 & time == 3
replace x1 = 16 if id == 3 & time == 4

replace x1 = 19 if id == 4 & time == 1
replace x1 = 20 if id == 4 & time == 2
replace x1 = 22 if id == 4 & time == 3
replace x1 = 23 if id == 4 & time == 4

replace x1 = 30 if id == 5 & time == 1
replace x1 = 31 if id == 5 & time == 2
replace x1 = 32 if id == 5 & time == 3
replace x1 = 33 if id == 5 & time == 4

replace x2 = 100 if id == 1 & time == 1
replace x2 = 102 if id == 1 & time == 2
replace x2 = 104 if id == 1 & time == 3
replace x2 = 106 if id == 1 & time == 4

replace x2 = 200 if id == 2 & time == 1
replace x2 = 201 if id == 2 & time == 2
replace x2 = 203 if id == 2 & time == 3
replace x2 = 205 if id == 2 & time == 4

replace x2 = 101 if id == 3 & time == 1
replace x2 = 103 if id == 3 & time == 2
replace x2 = 105 if id == 3 & time == 3
replace x2 = 107 if id == 3 & time == 4

replace x2 = 198 if id == 4 & time == 1
replace x2 = 200 if id == 4 & time == 2
replace x2 = 202 if id == 4 & time == 3
replace x2 = 204 if id == 4 & time == 4

replace x2 = 300 if id == 5 & time == 1
replace x2 = 301 if id == 5 & time == 2
replace x2 = 302 if id == 5 & time == 3
replace x2 = 303 if id == 5 & time == 4

order id time event x1 x2
sort id time
list, sepby(id)


ren x1 x 
ren x2 z 
mystacked_blad, rperiod(1) event(event) xvar(x z ) time(time) id(id)  seed(123)

*/
