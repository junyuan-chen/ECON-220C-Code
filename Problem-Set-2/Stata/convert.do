clear
capture log close
log using convert.log, replace
insheet using cigar.csv
save cigar.dta
log close
