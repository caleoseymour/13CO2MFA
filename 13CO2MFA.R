#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)

library('Rsolnp')
library('Rcpp')
library('RcppArmadillo')
library('truncdist')

## Place a file in the input fulder and then modify the fname variable to fit it
## (The framework around the model needs some work.)
#fname=NULL

dir.create('output/',showWarnings=FALSE)
dir.create('result/',showWarnings=FALSE)

Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
sourceCpp('model.cpp')

mcsim = 1000
nrestarts = 1
nsim = 100000

lims = data.frame(
r12 = c(0,100),
br1 = c(0,10),
r13 = c(0,100),
r14 = c(0,100))
rownames(lims) = c('min','max')

fin = args[1]
fout = args[2]

current_table = as.matrix(read.table(fin, header=TRUE, sep='\t',row.names=1))

inputs = apply(current_table, 2,
function(x)
{
    m = x[1]
    s = x[2]
    n = x[3]
    if (n > 1)
    {
        out = m + s * rtrunc(n=mcsim, spec='t', a=-Inf
        #qt(p=0.95,df=n-1)
        , b=Inf
        #qt(p=0.95,df=n-1)
        , df=n-1)
        out[out < 0] = 1e-6
        return(out)
    } else {
        #warning('Sample size is too small to estimate distribution.')
        return(replicate(mcsim, m))
    }
})
inputs = t(apply(inputs,1, function(x) return((x/sum(x))*sum(current_table[1,]))))
colnames(inputs) = colnames(current_table)
message('Running model.')
sink(fout)
    cat('r12\tbr1\tr13\tr14\terror\tCUE\tenergy\tobs[C1/CU]\tobs[C2/CU]\tobs[C3/CU]\tobs[C3/CU]\tobs[C5/CU]\tobs[C6/CU]\test[C1/CU]\test[C2/CU]\test[C3/CU]\test[C4/CU]\test[C5/CU]\test[C6/CU]\tr1\tr2\tr3\tr4\tr5\tr6\tr7\tr8\tr9\tr10\tr11\tr12\tr13\tr14\tr15\tbr1\tbr2\tbr3\tbr4\tbr5\tbr6\tbr7\tbr8\n')
sink()
for (i in 1:mcsim)
{
    setModelParam('input_values', matrix(inputs[i,],nrow=1))
    model_solution = gosolnp(fun=modelFunctionCpp
                             ,ineqfun=modelInequalityCpp
                             ,LB=unlist(lims[1,])
                             ,UB=unlist(lims[2,])
                             ,ineqLB=replicate(8, 0)
                             ,ineqUB=replicate(8,Inf)
                             ,n.sim=nsim
                             ,n.restarts=nrestarts
                             ,control=list(inner.iter=1000,outer.iter=2000,tol=1e-6))
    modelFunctionCpp(model_solution$pars)
    modelPrint(verbose=TRUE)
    sink(paste0(fout), append=TRUE)
        modelPrint(verbose=FALSE)
    sink()
}