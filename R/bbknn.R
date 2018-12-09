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
##'
##' @keywords graph network igraph mvtnorm simulation
##' @import reticulate
##' @importFrom stats prcomp
##' @export


bbknn <- function(data_matrix, batch, pca = TRUE, compute_pca = "python"){
    #import python modules with reticulate
    #reticulate::use_python("/usr/local/bin/python3")
    anndata <- reticulate::import("anndata",convert=FALSE)
    bbknn <- reticulate::import("bbknn", convert=FALSE)
    sc <- reticulate::import("scanpy.api",convert=FALSE)

    #set up annotation data for batches
    if(is.character(batch))  batch <- as.factor(batch)
    if(is.factor(batch))     batch <- as.numeric(batch)
    if(is.numeric(batch))    batch <- as.integer(batch)
    adata <- anndata$AnnData(X=data_matrix, obs=batch)
    
    #perform PCA
    if(pca){
        sc$tl$pca(adata)
        if(compute_pca == "python"){
            #use PCA computed in Python
            sc$pp$pca(data_matrix)
        } else {
            #use PCA computed in R
            adata$obsm$X_pca <- prcomp(data_matrix)$rotation
        }
        pca <- sc$pp$pca(data_matrix)
        adata$obsm$X_pca <- pca
    } else {
        #use full matrix
        reticulate::py_set_item(adata$obsm, name = "X_pca", value = reticulate::np_array(data_matrix))
    }
    
    
    #perform BBKNN to derive corrected components
    bbknn$bbknn(adata, batch_key=0)
    corrected_matrix <- as.matrix(adata$data)
    return(corrected_matrix)
}
