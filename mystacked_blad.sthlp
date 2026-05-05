{smcl}
{* *! version 1.0 26mar2026}{...}
{vieweralsosee "[R] help" "help help"}{...}
{viewerjumpto "Syntax" "mystacked_blad##syntax"}{...}
{viewerjumpto "Description" "mystacked_blad##description"}{...}
{viewerjumpto "Options" "mystacked_blad##options"}{...}
{viewerjumpto "Remarks" "mystacked_blad##remarks"}{...}
{viewerjumpto "Examples" "mystacked_blad##examples"}{...}

{title:Title}

{phang}
{bf:mystacked_blad} {hline 2} Build stacked matched samples for panel data with multiple treatment events

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:mystacked_blad}
[{it:if}]
[{it:in}],
{cmd:rperiod(}{it:integer}{cmd:)}
{cmd:event(}{it:varname numeric}{cmd:)}
{cmd:xvar(}{it:varlist}{cmd:)}
{cmd:time(}{it:varname numeric}{cmd:)}
{cmd:id(}{it:varname}{cmd:)}
[
{cmd:force}
{cmd:path_actual}
{cmd:by(}{it:varlist}{cmd:)}
{cmd:seed(}{it:#}{cmd:)}
{cmd:noyet(}{it:#}{cmd:)}
{cmd:allstack}
{cmd:onlynoyet}
]

{marker description}{...}
{title:Description}

{pstd}
{cmd:mystacked_blad} constructs stacked samples from panel data organized at the
{it:id} x {it:time} x {it:event} level. Units with {cmd:event>0} are treated,
while units with {cmd:event==0} are treated as never-treated controls.

{pstd}
By default, for each treated unit and event time, the command selects the nearest
eligible control using Mahalanobis distance based on the covariates supplied in
{cmd:xvar()}. These covariates are measured in the reference period
{cmd:event-rperiod}.

{pstd}
The command accommodates multiple treatment events because matching is performed
separately for each event time recorded in the data. If {cmd:by()} is specified,
eligible controls are restricted to the same group as the treated unit.

{pstd}
Alternatively, with {cmd:allstack}, the command returns the full stacked set of
treated-control comparisons instead of selecting a single nearest match.

{marker options}{...}
{title:Options}

{phang}
{cmd:rperiod(}{it:integer}{cmd:)} specifies how many periods before the event are
used to retrieve the matching covariates. For example, {cmd:rperiod(1)} uses
covariates measured one period before treatment.

{phang}
{cmd:event(}{it:varname numeric}{cmd:)} specifies the event-time variable.
Observations with {cmd:event>0} identify treated units and their event time.
Observations with {cmd:event==0} identify never-treated controls.

{phang}
{cmd:xvar(}{it:varlist}{cmd:)} specifies the covariates used to compute
Mahalanobis distance. These variables must be observed in the reference period
defined by {cmd:event-rperiod}.

{phang}
{cmd:time(}{it:varname numeric}{cmd:)} specifies the panel time variable.

{phang}
{cmd:id(}{it:varname}{cmd:)} specifies the panel unit identifier.

{phang}
{cmd:by(}{it:varlist}{cmd:)} restricts matching to occur within groups defined by
the variables in {cmd:by()}. For example, with {cmd:by(region sex)}, treated units
are matched only to controls in the same region and sex category.

{phang}
{cmd:noyet(}{it:#}{cmd:)} allows not-yet-treated units to be used as controls,
provided their own event occurs at least {it:#} periods after the treated unit's
event. Never-treated units remain eligible controls.

{phang}
{cmd:allstack} returns all eligible treated-control comparisons instead of
selecting the nearest control by Mahalanobis distance. This option is useful when
the researcher wants the full stacked comparison sample.

{phang}
{cmd:onlynoyet} is used with {cmd:allstack}. It keeps only comparisons involving
not-yet-treated controls and drops comparisons involving never-treated controls.

{phang}
{cmd:seed(}{it:#}{cmd:)} sets the random seed used to break ties when two or more
controls have the same Mahalanobis distance for a treated unit and event.

{phang}
{cmd:force} tells the command to proceed even if the expected number of
treated-control combinations exceeds 50,000,000. Without this option, execution
stops when the limit is exceeded.

{phang}
{cmd:path_actual} stores temporary files and the final matched dataset in the
current working directory. By default, the command uses Stata's temporary
directory.

{marker remarks}{...}
{title:Remarks}

{pstd}
The data must be uniquely identified by {cmd:id}, {cmd:time}, and {cmd:event}.
That is, there should be no duplicate observations within an
{it:id} x {it:time} x {it:event} cell.

{pstd}
By default, matching is one-to-one and is performed with replacement. A control
unit can therefore be matched to more than one treated unit or event.

{pstd}
When {cmd:allstack} is specified, the command does not compute nearest-neighbor
matches. Instead, it creates a stacked dataset containing all eligible treated and
control observations. In this case, the command also creates a {cmd:weight}
variable to account for repeated appearances of the same unit-event observation
across stacks.

{pstd}
The output dataset contains the original unit identifier, the original event time,
a new stack or pair identifier called {cmd:match_pair}, and a variable
{cmd:treated} indicating whether the observation corresponds to the treated or
control member of the stack.

{pstd}
When {cmd:noyet()} is not specified, controls are restricted to never-treated
units, defined by {cmd:event==0}. When {cmd:noyet()} is specified, not-yet-treated
units may also serve as controls if their treatment occurs sufficiently far in the
future.

{pstd}
If {cmd:seed()} is not specified, the command uses seed 123 for random
tie-breaking.

{marker examples}{...}
{title:Examples}

{pstd}
Nearest-neighbor matching using never-treated controls:

{phang2}
{cmd:. mystacked_blad, rperiod(1) event(event) xvar(x1 x2) time(time) id(id)}

{pstd}
Matching within groups:

{phang2}
{cmd:. mystacked_blad, rperiod(1) event(event) xvar(x1 x2) time(time) id(id) by(region sex)}

{pstd}
Allow not-yet-treated units to be used as controls if they are treated at least
two periods later:

{phang2}
{cmd:. mystacked_blad, rperiod(1) event(event) xvar(x1 x2) time(time) id(id) noyet(2)}

{pstd}
Return the full stacked set of eligible treated-control comparisons:

{phang2}
{cmd:. mystacked_blad, rperiod(1) event(event) xvar(x1 x2) time(time) id(id) allstack}

{pstd}
Return only stacked comparisons involving not-yet-treated controls:

{phang2}
{cmd:. mystacked_blad, rperiod(1) event(event) xvar(x1 x2) time(time) id(id) noyet(2) allstack onlynoyet}

{pstd}
Use the current working directory for temporary files:

{phang2}
{cmd:. mystacked_blad, rperiod(1) event(event) xvar(x1 x2) time(time) id(id) path_actual}

{pstd}
Set a seed for reproducible tie-breaking:

{phang2}
{cmd:. mystacked_blad, rperiod(1) event(event) xvar(x1 x2) time(time) id(id) seed(123)}

{pstd}
Proceed even if the number of possible treated-control combinations is large:

{phang2}
{cmd:. mystacked_blad if sample==1, rperiod(2) event(event) xvar(age income) time(year) id(firm_id) force}

{title:Author}

{pstd}
Bladimir Carrillo

{pstd}
Sao Paulo School of Economics-FGV

{pstd}
Email: {browse "mailto:blacaber8@hotmail.com":blacaber8@hotmail.com}

{psee}
Manual: {help help}
