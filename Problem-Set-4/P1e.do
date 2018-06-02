set obs 100
drawnorm x
gen y = x > 0
logit y x
