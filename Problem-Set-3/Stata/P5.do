clear
local beta -0.5
cap postclose temp
postfile temp beta_hat_2sls beta_hat_iv se_hat_2sls se_hat_iv ///
           using iv_estimate,replace
    /* declares the variable names and the filename of a (new) Stata dataset
where
    results are to be stored.  */
timer clear 1
timer on 1

forvalues i = 1(1)1000 { /* perform the experiment 1000 times */
    drop _all /* drop all the variables in memory */
    quietly set obs 1000
    gen z = rnormal(0,1)
    gen xc = ln(2)*z
    gen xd = rnormal(0,1)
    gen u = 0.5*xd + rnormal(0,1)
    gen x = xc + xd
    gen y = 1 -`beta'*x + u /* try -0.5,0.5 and 1 */
    qui reg x z
    qui predict x_hat, xb
    qui reg y x_hat
    scalar beta_hat_2sls = _b[x_hat]
    scalar se_hat_2sls = _se[x_hat]
    qui ivreg y (x=z)
    scalar beta_hat_iv = _b[x]
    scalar se_hat_iv = _se[x]
    post temp (beta_hat_2sls) (beta_hat_iv) (se_hat_2sls) (se_hat_iv)
}
timer off 1
postclose temp
use iv_estimate, clear
sum
timer list 1
