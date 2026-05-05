# mystacked_blad

`mystacked_blad` is a Stata command that builds a stacked matched sample using Mahalanobis distance.

The command is designed for panel data organized at the id × time × event level. For each treated unit and event time, it constructs eligible control sets—optionally within user-defined groups—and either selects the nearest control using user-supplied covariates measured in a pre-treatment reference period or, with allstack, returns the full stacked set of treated-control comparisons. Controls can be restricted to never-treated units or expanded to include not-yet-treated units through the noyet option, with optional filtering via onlynoyet. The command naturally accommodates settings with multiple treatment events, as matching is performed separately for each event time recorded in the data.

## Installation

Install directly from GitHub in Stata with:

```stata
net install mystacked_blad, from("https://raw.githubusercontent.com/blacaber8/mystacked_blad/main/") replace
