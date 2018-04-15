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

* Question 1

*reg log_vio shall, r
/*
esttab using P5Q1.tex, label nostar replace booktabs
reg log_mur shall, r
esttab using P5Q1.tex, label nostar append booktabs
reg log_rob shall, r
esttab using P5Q1.tex, label nostar append booktabs
*/
eststo clear
eststo: reg log_vio shall, r
eststo: reg log_mur shall, r
eststo: reg log_rob shall, r
esttab using P5Q1.tex, se r2 nostar replace booktabs
