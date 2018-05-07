clear
capture log close
log using convert.log, replace
insheet using cigar.csv
rename logconsumption logC
rename lnprice logP
rename logincome logY
rename logminimumpriceatneighboringstat logPn
save cigar.dta, replace
log close
