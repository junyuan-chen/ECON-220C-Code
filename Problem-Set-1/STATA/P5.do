clear
*set more off;
*set matsize 300;
capture log close
/* ps3 example using handguns data */
log using shall.log, replace
use handguns.dta
desc
summarize
gen log_vio=log(vio)
gen log_mur=log(mur)
gen log_rob=log(rob)

* Part I

eststo clear
eststo: reg log_vio shall, r
eststo: reg log_mur shall, r
eststo: reg log_rob shall, r
esttab using P5I.tex, se r2 nostar replace booktabs

* Part II

eststo clear
eststo: reg log_vio shall incarc_rate density pop pm1029 avginc, r
eststo: reg log_mur shall incarc_rate density pop pm1029 avginc, r
eststo: reg log_rob shall incarc_rate density pop pm1029 avginc, r
esttab using P5II.tex, se r2 nostar replace booktabs

* Part IV

tab state, gen(statedummy)
/* column 1 in the table */
reg log_vio shall incarc_rate density pop pm1029 avginc, r
/* column 2 in the table */
reg log_vio shall incarc_rate density pop pm1029 avginc statedummy*, r
testparm statedummy*
/* if you want to compute standard error that is robust to the time series
correlation in uit, you can u se the following commands. */
reg log_vio shall incarc_rate density pop pm1029 avginc statedummy*,cluster(state) r
/* you can also use xtreg here */
*Dep=ln(violence)
