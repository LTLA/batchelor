# Adding a test to flush out any uncontrolled parallelization.
library(BiocParallel)
failgen <- setRefClass("FailParam", 
    contains="BiocParallelParam",     
    fields=list(),
    methods=list())

FAIL <- failgen()
register(FAIL) 

library(DelayedArray)
setAutoBPPARAM(FAIL)
