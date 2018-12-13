##' @name bbknn State Matrix
##' @rdname bbknn
##'
##' @title Run bbknn clustering algorithm
##'
##' @description Implements the bbknn clustering algorithm in R using reticulate to run the Python version. Requires the python "bbknn" and "igraph" modules to be installed. Returns a vector of partition indices.
##'
##' @param data_matrix A matrix (genes x samples or cells) for expression data
##' @param batch An integer vector of batches to correct for (converts factors or numeric vectors)
##' @param pca whether to compute pca (defaults to TRUE) or apply correction to the raw matrix (FALSE) 
##' @param compute_pca whether to compute PCA in Python (defaults to TRUE, requires scanpy library) or with R functions (FALSE)
##' @param nPcs number of principal components to compute (defaults to 50 if more than 50 genes)
##'
##' @return returns a list with the following components \item{corrected matrix}{matrix of data corrected by the BBKNN (batch based K nearest neighbours)}\item{pca}{principal components(matrix with row for every sample and column for each component)}\item{tsne}{t-distributed stochastic neighbour embedding (matrix with row for every sample)}\item{umap}{uniform manifold approximation and projection (matrix with row for every sample)}
##'
##' @keywords graph network igraph mvtnorm simulation
##' @import reticulate
##' @importFrom stats prcomp
##' @export


bbknn <- function(data_matrix, batch, pca = TRUE, compute_pca = "python", nPcs = NULL){
    #import python modules with reticulate
    if(!is.matrix(data_matrix)){
        warning("matrix expected for data_matrix")
        data_matrix <- as.matrix(data_matrix)
    } 
    if(is.null(nPcs)) nPcs <- min(50, nrow(data_matrix), ncol(data_matrix))
    if(nPcs > nrow(data_matrix)){
        warning("number of genes less than nPcs")
        print(paste("using", nrow(data_matrix), "components"))
        nPcs <- nrow(data_matrix)
    }
    #reticulate::use_python("/usr/local/bin/python3")
    anndata <- reticulate::import("anndata",convert=FALSE)
    bbknn <- reticulate::import("bbknn", convert=FALSE)
    sc <- reticulate::import("scanpy.api",convert=FALSE)

    #set up annotation data for batches
    if(is.character(batch))  batch <- as.factor(batch)
    if(is.factor(batch))     batch <- as.numeric(batch)
    if(is.numeric(batch))    batch <- as.integer(batch)
    
    #perform PCA
    if(pca){
        #sc$tl$pca(adata)
        if(compute_pca == "python"){
            #use PCA computed in Python
            pca <- sc$pp$pca(t(data_matrix))
        }else if(compute_pca != "python"){
            #use PCA computed in R
            print("test")
            pca <- reticulate::r_to_py(t(prcomp(data_matrix)$x[1:nPcs,]))
        }
        adata <- anndata$AnnData(X=pca, obs=batch)
        sc$tl$pca(adata)
        adata$obsm$X_pca <- pca
    } else {
        #use full matrix
        adata <- anndata$AnnData(X=t(data_matrix), obs=batch)
        sc$tl$pca(adata)
    }
    #perform BBKNN to derive corrected components
    bbknn$bbknn(adata, batch_key = 0)
    corrected_matrix <- t(reticulate::py_to_r(adata$X))
    sc$tl$pca(adata)
    pca <- reticulate::py_to_r(adata$obsm$X_pca)
    sc$tl$tsne(adata)
    tsne <- reticulate::py_to_r(adata$obsm$X_tsne)
    sc$tl$umap(adata)
    umap <- reticulate::py_to_r(adata$obsm$X_umap)
    output <- list(corrected_matrix = corrected_matrix, pca = pca, tsne = tsne, umap = umap)
    return(output)
}




