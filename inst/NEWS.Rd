\name{batchelornews}
\title{batchelor News}
\encoding{UTF-8}

\section{Version 1.4.0}{\itemize{
\item Support the use of arbitrary design matrices in \code{regressBatches()}.

\item Allow lists of objects to be directly passed into the \code{...} for many functions.

\item Added the \code{clusterMNN()} function for performing MNN on cluster centroids.

\item Added \code{get.variance=} option to \code{fastMNN()} to return variance explained from PCA.
Support \code{d=NA} to skip the PCA step altogether.
}}

\section{Version 1.2.0}{\itemize{
\item Deprecated rotate.all= in favour of get.all.genes= in multiBatchPCA().

\item Switched BSPARAM= to use IrlbaParam(deferred=TRUE) by default in fastMNN(),
so that the default behaviour is actually fast.

\item Deprecated auto.order= in favor of merge.order= and auto.merge= in fastMNN() and mnnCorrect().
Automatic merging now detects potential tree-based merges.
Merge trees can also be specified as input.

\item Added the correctExperiments() function to cbind the original assays alongside the merged values.

\item Added the subset.row= argument to cosineNorm() for in-place subsetting.

\item Added batch= and preserve.single= arguments to multiBatchNorm().
Standardized behavior of subset.row= by adding a normalize.all= argument.

\item Added the regressBatches() function for correction via standard linear regression. 

\item Added the prop.k= argument in all MNN-related functions,
to allow the value of k to adapt asymmetrically to the size of each batch.
}}

\section{Version 1.0.0}{\itemize{
\item
New package batchelor, for batch correction of single-cell (RNA sequencing) data.
}}
