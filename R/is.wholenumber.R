# -------------------------------------
# RCropStat Utilities: Function
# -------------------------------------
# isWholenumber/is.wholenumber:
# Reference: FROM FULLREFMAN.PDF
# 
# Author: Alaine Gulles
# -------------------------------------

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) { 
	abs(x - round(x)) < tol 
}
