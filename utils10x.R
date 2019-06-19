addstringToColnames = function(object, name=NULL){
  if (is.null(name)) name = deparse(substitute(object))
  colnames(object@mat) = paste0(name,"_",colnames(object@mat)) 
  return(object)
}

obtainfilteredCountMat = function(object, min.trans = 100, collapsetogenes=TRUE){
  #object = cellranger object
  Count.mat = object@mat
  
  # filter non-expressed genes and cells with fewer than min.trans transcripts
  genes.use = apply(Count.mat,1,sum) > 0 # Expressed genes
  Count.mat = Count.mat[genes.use,]
  Count.mat = Count.mat[,apply(Count.mat,2,sum) > min.trans]
  
  if (collapsetogenes){
     which.dup = names(table(object@gene_symbols[genes.use]))[table(object@gene_symbols[genes.use])>1]
     temp.mat = data.frame()
     geneNames = object@gene_symbols[genes.use]
     for (gene in which.dup){
       rows1 = which(geneNames %in% gene)
       geneNames = geneNames[-rows1]
       temp.mat = Matrix::colSums(Count.mat[rows1,])
       Count.mat = Count.mat[-rows1,]
       Count.mat = rbind(Count.mat, temp.mat)
       rownames(Count.mat)[nrow(Count.mat)] = gene
       geneNames = c(geneNames, gene)
     }
  }
  

  rownames(Count.mat) = geneNames
  return(Count.mat)
}