clear
capture postclose tempid
postfile tempid beta1 beta2 beta3 se1 se2 se3 using mydata.dta, replace
forvalues i = 1(1)100 {
drop _all
quietly set obs 50000
/************* DGP ********************/
gen e1 = rnormal()
gen e2 = rnormal()
*gen z = rnormal()
gen z = 5*rnormal()
gen x1 = z + e1
gen x2 = rnormal()
gen x = x1 + x2
gen u = x2 + e2
gen y_latent = x + u
gen d = 0
qui replace d = 1 if y_latent > 0
gen y = y_latent*d
/********** the first probit regression ************/
quietly probit d x,r
scalar beta1 = _b[x]
scalar se1 = _se[x]
/********** the second probit regression ***********/
quietly reg x z
predict e, resid
quietly probit d x e
scalar beta2 = _b[x]
scalar se2 = _se[x]
/********** the third probit regression ***********/
quietly reg x z
predict x_hat, xb
quietly probit d x_hat
scalar beta3 = _b[x]
scalar se3 = _se[x]
post tempid (beta1) (beta2) (beta3) (se1) (se2) (se3)
}
postclose tempid
use mydata.dta, clear
sum
