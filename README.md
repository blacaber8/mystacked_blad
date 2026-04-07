# mystacked_blad

`mystacked_blad` is a Stata command that builds a stacked matched sample using Mahalanobis distance.

The command is designed for panel data organized at the `id × time x event` level. For each treated unit, it finds the nearest control based on user-supplied covariates measured in a reference period before treatment. The command naturally accommodates settings with multiple treatment events, as matching is performed separately for each event time recorded in the data.

## Installation

Install directly from GitHub in Stata with:

```stata
net install mystacked_blad, from("https://raw.githubusercontent.com/blacaber8/mystacked_blad/main/") replace
