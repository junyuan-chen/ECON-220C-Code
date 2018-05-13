* ECON 220C Problem Set 2
* Problem 4 (a)-(e)
* Junyuan Chen

clear
eststo clear
capture log close
log using P4.log, replace
use cigar.dta
desc
xtset state year, yearly

* Part (a)
eststo: reg logC L.logC logP logY logPn, r

* Part (b)
eststo: reg logC L.logC logP logY logPn i.year, r

* Part (c)
eststo: xtreg logC L.logC logP logY logPn, r fe

* Part (d)
eststo: xtreg logC L.logC logP logY logPn i.year, r fe
testparm i.year

* Part (e)
eststo: xi: xtivreg logC (L.logC = L2.logC) logP logY logPn i.year, fd
esttab using Tables/P4a-e.tex, se nostar drop (*year* _cons) booktabs replace

* Part (e) CORRECT WAY
xi: ivregress gmm D.logC (LD.logC = L2.logC) D.logP D.logY D.logPn i.year

log close
