{smcl}
{* *! version 1.0 26mar2026}{...}
{vieweralsosee "[R] help" "help help"}{...}
{viewerjumpto "Syntax" "mystacked_blad##syntax"}{...}
{viewerjumpto "Description" "mystacked_blad##description"}{...}
{viewerjumpto "Options" "mystacked_blad##options"}{...}
{viewerjumpto "Examples" "mystacked_blad##examples"}{...}

{title:Title}

{phang}
{bf:mystacked_blad} {hline 2} Build a stacked matched sample using Mahalanobis distance

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:mystacked_blad}
[{it:if}]
[{it:in}],
{cmd:rperiod(}{it:integer}{cmd:)}
{cmd:event(}{it:varname}{cmd:)}
{cmd:xvar(}{it:varlist}{cmd:)}
{cmd:time(}{it:varname}{cmd:)}
{cmd:id(}{it:varname}{cmd:)}
[{cmd:force}
{cmd:path_actual}
{cmd:by(}{it:varlist}{cmd:)}
{cmd:seed(}{it:#}{cmd:)}]

{marker description}{...}
{title:Description}

{pstd}
{cmd:mystacked_blad} constructs a stacked matched sample from panel data organized at the {it:id} x {it:time} x {it:event} level.
For each treated unit, the command finds the nearest control using Mahalanobis distance computed from the variables supplied in {cmd:xvar()} and measured in period {cmd:event-rperiod}.

{pstd}
If {cmd:by()} is specified, matching is restricted to occur within the groups defined by the variables in {cmd:by()}.
For example, with {cmd:by(region sex)}, a treated unit is matched only to controls in the same region and sex category.

{marker options}{...}
{title:Options}

{phang}
{cmd:rperiod(}{it:integer}{cmd:)} specifies how many periods before the event should be used to retrieve matching covariates.

{phang}
{cmd:event(}{it:varname}{cmd:)} specifies the event-time variable. Units with {cmd:event>0} are treated; units with {cmd:event==0} are controls.

{phang}
{cmd:xvar(}{it:varlist}{cmd:)} specifies the covariates used in the Mahalanobis distance calculation.

{phang}
{cmd:time(}{it:varname}{cmd:)} specifies the panel time variable.

{phang}
{cmd:id(}{it:varname}{cmd:)} specifies the panel unit identifier.

{phang}
{cmd:force} tells the command to proceed even if the expected number of treated-control combinations exceeds 50,000,000.

{phang}
{cmd:path_actual} tells the command to use the current working directory to store temporary files and the final matched dataset.

{phang}
{cmd:by(}{it:varlist}{cmd:)} restricts matching to occur only within groups defined by the variables in {cmd:by()}.

{phang}
{cmd:seed(}{it:#}{cmd:)} sets the random seed used to break ties when two or more control units have the same Mahalanobis distance for a given treated unit and event.

{title:Remarks}

{pstd}
The data must be uniquely identified by {cmd:id}, {cmd:time} and {cmd:event}. In other words, the dataset must be at the {it:id} x {it:time} x {it:event} level.

{pstd}
Matching is one-to-one and is performed with replacement.

{pstd}
If {cmd:by()} is specified, treated units are matched only to controls within the same {cmd:by()} group.

{pstd}
If {cmd:path_actual} is specified, files are written to the current working directory instead of Stata's temporary directory.

{pstd}
If {cmd:seed()} is specified, tie-breaking becomes reproducible across runs.

{marker examples}{...}
{title:Examples}

{phang2}{cmd:. mystacked_blad, rperiod(1) event(event) xvar(x1 x2) time(time) id(id)}

{phang2}{cmd:. mystacked_blad, rperiod(1) event(event) xvar(x1 x2) time(time) id(id) by(region sex)}

{phang2}{cmd:. mystacked_blad, rperiod(1) event(event) xvar(x1 x2) time(time) id(id) path_actual}

{phang2}{cmd:. mystacked_blad, rperiod(1) event(event) xvar(x1 x2) time(time) id(id) seed(123)}

{phang2}{cmd:. mystacked_blad if sample==1, rperiod(2) event(event) xvar(age income) time(year) id(firm_id) force}

{title:Author}

{pstd}
Bladimir Carrillo

{pstd}
Sao Paulo School of Economics-FGV

{pstd}
Email: {browse "mailto:blacaber8@hotmail.com":blacaber8@hotmail.com}

{psee}
Manual: {help help}
