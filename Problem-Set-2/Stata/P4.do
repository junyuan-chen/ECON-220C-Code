* ECON 220C Problem Set 2
* Problem 4
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
esttab using Tables/P4a.tex, se nostar booktabs replace

* Part (b)
eststo: reg logC L.logC logP logY logPn i.year, r
esttab using Tables/P4b.tex, se nostar drop (*year*) booktabs replace

* Part (c)
eststo: xtreg logC L.logC logP logY logPn, r fe
esttab using Tables/P4c.tex, se nostar booktabs replace

* Part (d)
eststo: xtreg logC L.logC logP logY logPn i.year, r fe
esttab using Tables/P4d.tex, se nostar drop (*year*) booktabs replace
testparm i.year

* Part (e)
eststo: xi: xtivreg logC (L.logC = L2.logC) logP logY logPn i.year, fd
esttab using Tables/P4e.tex, se nostar drop (*year*) booktabs replace

* Part (g)
