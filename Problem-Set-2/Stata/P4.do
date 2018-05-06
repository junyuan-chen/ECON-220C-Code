* ECON 220C Problem Set 2
* Problem 4
* Junyuan Chen

clear
eststo clear
capture log close
log using P4.log, replace
use cigar.dta
desc


* Part (a)
gen logC_lag = logC[_n-1]
eststo: reg logC logC_lag logP logY logPn, r
esttab using Tables/P4a.tex, se nostar booktabs replace

* Part (b)
eststo: reg logC logC_lag logP logY logPn i.year, r
esttab using Tables/P4b.tex, se nostar drop (*year*) booktabs replace

* Part (c)
xtset state
eststo: xtreg logC logC_lag logP logY logPn, r fe
esttab using Tables/P4c.tex, se nostar booktabs replace

* Part (d)
eststo: ivregress logC logC_lag logP logY logPn i.year, r fe
esttab using Tables/P4d.tex, se nostar drop (*year*) booktabs replace
