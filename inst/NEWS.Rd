\name{batchelornews}
\title{batchelor News}
\encoding{UTF-8}

\section{Version 1.2.0}{\itemize{
\item Deprecated rotate.all= in favour of get.all.genes= in multiBatchPCA().

\item Switched BSPARAM= to use IrlbaParam(deferred=TRUE) by default in fastMNN(),
so that the default behaviour is actually fast.

\item Deprecated auto.order= in favor of merge.order= and auto.merge= in fastMNN() and mnnCorrect().
Automatic merging now detects potential tree-based merges.
Merge trees can also be specified as input.

\item Added the correctExperiments() function to cbind the original assays alongside the merged values.

\item Added the subset.row= argument to cosineNorm() for in-place subsetting.
}}

\section{Version 1.0.0}{\itemize{
\item
New package batchelor, for batch correction of single-cell (RNA sequencing) data.
}}
