# re-built the DWD packages based on the source code
## introduction
- the DWD packages is created with [Hanwen,H ,et.](https://academic.oup.com/bioinformatics/article/28/8/1182/195830?login=false) and for providing the implementation of distance weighted discrimination (i.e. DWD) using an interior point method for the solution of second order cone programming problems
- the problem is DWD is never updated for R 3.1.0, and there aren't any R packages which can perform the DWD to get weighted matrix, so it's a good idea to modify the source code to fit the new R version like R 4.2.0
## preparation
- DWD source packages [DWD](https://cran.r-project.org/src/contrib/Archive/DWD/)
## main problems
### kdwd.R
- ybuf is not found
  - caused by the if-condition that require input [y] as factor, it will skip if not input the factor, which will cause not assignment to yubf
### checkdepconstr.R
- different row number between UU and tmp when col-binding
  - caused by the tmp matrix change to vector automatically when the col=1
### qprod.R
- spMatrix parameter are wrong of x
  - caused by x here are matrix but not numeric, so need to change the class of x to numeric
### other R
- cBind and rBind are wrong function
  - may be the cBind and rBind are replaced with cbind and rbind after R updated

## the new package is in the champeil/play_software/R_package/DWD/DWD_0.11.tar.gz
