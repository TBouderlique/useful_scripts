#Replace the original function that was somehow removed from my package

p2.generate.go.web  <- function (gene.names, egALIAS2EG = NULL, egGO2ALLEGS = NULL, n.cores = 1) {
  if (!requireNamespace("GO.db", quietly = TRUE)) {
    stop("Package \"GO.db\" needed for this function to work. Please install it with `BiocManager::install('GO.db')`.", call. = FALSE)
  }
  
  if (is.null(egALIAS2EG)) {
    stop("egALIAS2EG cannot be null, it has to be an object like org.Hs.egALIAS2EG")
  }
  
  if (!is.character(gene.names)) {
    stop("gene.names needs to be a character vector of gene names")
  }
  
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
    stop("Package \"AnnotationDbi\" needed for this function to work. Please install it with `BiocManager::install('AnnotationDbi')`.", call. = FALSE)
  }
  
  ids <- unlist(mclapply(AnnotationDbi::mget(gene.names, egALIAS2EG, ifnotfound = NA), function(x) x[1], mc.cores = n.cores))
  rids <- names(ids)
  names(rids) <- ids
  go.env <- AnnotationDbi::eapply(egGO2ALLEGS, function(x) as.character(na.omit(rids[x])))
  go.env <- go.env[unlist(lapply(go.env, length)) > 5]
  
  ## Put the GO Term annotation generated in a format suitable for the web interface
  nms <- names(go.env)
  names(nms) <- nms
  geneSets <- lapply(nms, function(x) {
    list(
      properties = list(
        locked = TRUE,
        genesetname = x,
        shortdescription = GO.db::GOTERM[[x]]@Term
      ),
      genes = c(go.env[[x]])
    )
  })
  
  invisible(geneSets)
}


###################################################################### Differential expression between 2 conditions in each cluster

#x-> your seurat object
#clus.ter-> first metadata you wanna use in your comparison
#condition-> second metadata you wanna combine in your analysis
#ncore= number of cores you wanna use in the analysis


FindMarker.cluster.condition <- function(x,clus.ter, condition,ncore=5){
  library(Seurat)
  library(BiocParallel)
  library(dplyr)
  ##create metadata for comparison
  
  Idents(object=x) <- clus.ter
  x@meta.data$clus.ter<-x@active.ident
  
  Idents(object=x) <- condition
  x@meta.data$condition<-x@active.ident
  
  x@meta.data$comparison<-paste0(x@meta.data$clus.ter,'_',x@meta.data$condition)
  
  Idents(object=x) <- "comparison"
  
  condition1= unique(x@meta.data$condition)[1]
  condition2= unique(x@meta.data$condition)[2]
  
  M=length(unique(x@meta.data$clus.ter))
  N=M-1
  
  ##base function for comparison
  FindMarker.wrapper<- function(y){
    list<-matrix()
    FindMarkers(x,ident.1=paste0(y,'_',condition1),ident.2 =paste0(y,'_',condition2), only.pos = TRUE, min.pct=0.1)
  }
  #start actual comparison
  Markers <- bplapply(0:N, FindMarker.wrapper,BPPARAM=MulticoreParam(ncore))
  
  for (i in unique(x@meta.data$clus.ter)){
    Markers[[as.numeric(i)+1]]$cluster<-i
  }
  #make a simple table of all the comparisons
  markers <- dplyr::bind_rows(Markers)
  
  
  gene_list<-matrix()
  genes<-strsplit(rownames(markers),'[.]')
  for (i in 1:length(genes)){
    gene_list<-append(gene_list, genes[[i]][1])
  }
  
  gene_list<-gene_list[-1]
  markers$gene<-gene_list
  
  
  top5 <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
  top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  top20<- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
  
  list_diff_1<-list(markers, top5, top10, top20)
  names(list_diff_1)<- c(paste0('gene_list_',condition1,'_VS_',condition2,'_in_',clus.ter), 
                         paste0('top5_',condition1,'_VS_',condition2,'_in_',clus.ter), 
                         paste0('top10_',condition1,'_VS_',condition2,'_in_',clus.ter), 
                         paste0('top20_',condition1,'_VS_',condition2,'_in_',clus.ter))
  
  #return(list_diff_1)
  
  ####reverse
  condition1= unique(x@meta.data$condition)[2]
  condition2= unique(x@meta.data$condition)[1]
  
  ##base function for comparison
  FindMarker.wrapper<- function(y){
    list<-matrix()
    FindMarkers(x,ident.1=paste0(y,'_',condition1),ident.2 =paste0(y,'_',condition2), only.pos = TRUE, min.pct=0.1)
  }
  #start actual comparison
  Markers <- bplapply(0:N, FindMarker.wrapper,BPPARAM=MulticoreParam(ncore))
  
  for (i in unique(x@meta.data$clus.ter)){
    Markers[[as.numeric(i)+1]]$cluster<-i
  }
  #make a simple table of all the comparisons
  markers <- dplyr::bind_rows(Markers)
  
  
  gene_list<-matrix()
  genes<-strsplit(rownames(markers),'[.]')
  for (i in 1:length(genes)){
    gene_list<-append(gene_list, genes[[i]][1])
  }
  
  gene_list<-gene_list[-1]
  markers$gene<-gene_list
  
  
  top5 <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
  top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  top20<- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
  
  list_diff_2<-list(markers, top5, top10, top20)
  names(list_diff_2)<- c(paste0('gene_list_',condition1,'_VS_',condition2,'_in_',clus.ter), 
                         paste0('top5_',condition1,'_VS_',condition2,'_in_',clus.ter), 
                         paste0('top10_',condition1,'_VS_',condition2,'_in_',clus.ter), 
                         paste0('top20_',condition1,'_VS_',condition2,'_in_',clus.ter))
  
  res<-list()
  res$list1<-list_diff_1
  res$list2<-list_diff_2
  
  names(res)<- c(paste0(condition2,'_VS_',condition1), 
                paste0(condition1,'_VS_',condition2))
                         
  
  return(res)
  
}



####################################### Find markers per type fo metadata

Find.markers.perCluster<- function(x,project='',cluster,ncore=5){
  
  library(Seurat)
  library(BiocParallel)
  library(dplyr)
  
  Idents(x) <- cluster
  M=unique(x@active.ident)

 FindMarker.wrapper <- function(y){
    list<-matrix()
    FindMarkers(x,ident.1=y, only.pos = TRUE, min.pct=0.1)
  }
  
  Markers <- bplapply(M, FindMarker.wrapper,BPPARAM=MulticoreParam(ncore))
    
   for (i in 1:length(unique(x@active.ident))){
    Markers[[as.numeric(i)]]$cluster<-unique(x@active.ident)[i]
  }
  markers <- dplyr::bind_rows(Markers)  
  
  gene_list<-matrix()
  genes<-gsub(pattern = "\\.\\..*", replacement = "",x = rownames(markers))
  for (i in 1:length(genes)){
    gene_list<-append(gene_list, genes[[i]][1])
  }
  gene_list<-gene_list[-1]
  markers$gene<-gene_list
    
  
  
  top5 <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
  top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  top20<- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
  
  list_diff_cluster<-list(markers, top5, top10, top20)
  names(list_diff_cluster)<- c(paste0(project,'diff_allgenes'),
                               paste0(project,'diff_top5'),
                               paste0(project,'diff_top10'),
                               paste0(project,'diff_top20')
                              )
  
  
  return(list_diff_cluster)
  
}


############################################## human to mouse genes

convertHumanToMouseGeneList<- function(x){
  library("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  mousex <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(mousex))
  return(mousex)
}

###############################################################mouse to human genes

convertMouseToHumanGeneList<- function(x){
  library("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = mouse, attributesL = c("mgi_symbol"), martL = human, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}


#############################################################Ensembl to mouse genes
EnsemblToMouseGene<-function(x){
  
  s<-x
  library('biomaRt')
ensembl <- useEnsembl(biomart = "ensembl")
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
genes <- rownames(s@assays$RNA@counts)
ensembl <- useDataset(dataset = "mmusculus_gene_ensembl", mart = ensembl)
x<-getBM(attributes = c('external_gene_name', 'ensembl_gene_id'),    # info to be recovered
         filters = 'ensembl_gene_id',                                #restrict output to ensembl gene ids
         values = rownames(s@assays$RNA@counts),                     #ensembl values to recover
         mart = ensembl,
         useCache = F)
#replace values in rownames
require(plyr)

rownames(s@assays$RNA@counts)<-mapvalues(rownames(s@assays$RNA@counts), 
                                         from=x$ensembl_gene_id, 
                                         to=x$external_gene_name)

rownames(s@assays$RNA@data)<-mapvalues(rownames(s@assays$RNA@data), 
                                       from=x$ensembl_gene_id, 
                                       to=x$external_gene_name)

rownames(s@assays$RNA@data)<-make.unique(rownames(s@assays$RNA@data))
rownames(s@assays$RNA@counts)<-make.unique(rownames(s@assays$RNA@counts))

s[["RNA"]]@meta.features <- data.frame(row.names = rownames(s[["RNA"]]))

return(s)
}

################################################Make seurat basic
MakeSeurat.wrapper<-function(x, sample){
  rownames(x)<-make.unique(rownames(x))
  
  seurat<-CreateSeuratObject(x,min.cells = 2, min.features = 200)
  
  
  mito_genes <- grep(pattern = "^mt-", x = rownames(x = seurat@assays$RNA@data), value = TRUE,ignore.case = T)
  percent_mito <- Matrix::colSums(seurat@assays$RNA@data[mito_genes, ])/Matrix::colSums(seurat@assays$RNA@data)
  
  seurat@meta.data$percent_mito=percent_mito
  
  q1<-VlnPlot( seurat, features = c("nCount_RNA", "nFeature_RNA", "percent_mito"),pt.size = 0.5)
  
  base<-getwd()
  dir.create(paste0(base,'/exports/QC/'),showWarnings = F)
  ggsave(paste0(base,'/exports/QC/',sample,'.png'),cowplot::plot_grid(q1,nrow = 1), width = 15,height = 15, units = 'cm')
  
  
  
  return(seurat)
}



################################# Find doublets



get.scrublet.scores <- function(mat, min.molecules.per.gene=10) {
  require(data.table)
   # write out a CSV file
  tf <- tempfile()
  dtf <- paste(tf,'doubletScores',sep='.')
  dtp <- paste(tf,'predicted_doublets',sep='.')
  dt <- data.table(as.matrix(t(mat[Matrix::rowSums(mat)>=min.molecules.per.gene,])))
  data.table::fwrite(dt,file=tf)
  cmd <- paste("/home/tbou/miniconda3/envs/RNAvelo/bin/python3.8 -c 'import sys; import pandas; import scrublet; df = pandas.read_csv(\"",tf,"\"); scrub = scrublet.Scrublet(df); doublet_scores, predicted_doublets = scrub.scrub_doublets(); pandas.DataFrame(doublet_scores).to_csv(\"",dtf,"\"); pandas.DataFrame(predicted_doublets).to_csv(\"",dtp,"\");'",sep='')
  
  ### Be aware of which python you have installed scrublet at! Also, you can modify the function here to save the histogram of scrublet scores.
  
  print(cmd)
  tmp <- system(cmd, intern=T)
  #system(cmd);
  x <- as.data.frame(data.table::fread(dtf,sep=',',header=F,skip=1))
  y <- as.data.frame(data.table::fread(dtp,sep=',',header=F,skip=1))
  x <- as.numeric(as.data.frame(data.table::fread(dtf,sep=',',header=F,skip=1))[,2])
  y <- as.numeric(as.data.frame(data.table::fread(dtp,sep=',',header=F,skip=1))[,2])
  names(x)<-colnames(mat)
  names(y)<-colnames(mat)
  file.remove(tf)
  file.remove(dtf)
  file.remove(dtp)
  
  
  z<-rbind(x,y)
  
  return(z)
}



doublets_col <- function(sample) {    ###  sample is Seurat object
  mat <- seurat@assays$RNA@counts
  # get scores from scrublet
  doublet_scores <- get.scrublet.scores(mat)
  sample@meta.data$doublet_score <- doublet_scores[1,]
  sample@meta.data$predicted_doublet <- as.logical(doublet_scores[2,])
  return(sample)
}


######################for Conos integration
quick.seurat.From.list<-function(x) {
  for (i in 1:length(x)) {
  x[[i]] <- NormalizeData(x[[i]], verbose = FALSE)
  x[[i]] <- FindVariableFeatures(x[[i]], selection.method = "vst", 
                                           nfeatures = 4000, verbose = FALSE)
  
  x[[i]] <- ScaleData(x[[i]], verbose = FALSE)
  x[[i]] <- RunPCA(x[[i]], npcs = 30, verbose = FALSE)
  x[[i]] <- RunUMAP( x[[i]], reduction = "pca", dims = 1:30)
  #x[[i]] <- RunTSNE( x[[i]], reduction = "pca", dims = 1:30)
  
  
  }
  return(x)
 }



Conos_integration<-function(x,n.cores=5, k=30,k.self=5, space='PCA', ncomps=30, n.odgenes=2000,matching.method='mNN',metric='angular',score.component.variance=F){
  require(conos)
  
  con <- Conos$new(x, n.cores=n.cores)
  con$buildGraph(k=k, k.self=k.self, space=space, ncomps=ncomps, 
                 n.odgenes=n.odgenes, matching.method=matching.method, 
                 metric=metric, score.component.variance=score.component.variance, 
                 verbose=TRUE)
  
  con$findCommunities(method=leiden.community, resolution=1)
  
  con$embedGraph(method='UMAP')
  conos_umap<-con$embedding
  colnames(conos_umap)<-c('UMAP_conos_1','UMAP_conos_1')
  
  return(conos_umap)
  
}


 
                     
                                  
                                  
###############################################################Cytotrace#####################
                                  
  CytoTRACE <- function(mat, batch = NULL, enableFast = TRUE,
                      ncores = 1,subsamplesize = 1000){

  range01 <- function(x){(x-min(x))/(max(x)-min(x))}

  #inputs
  a1 <- mat
  a2 <- batch
  if(ncol(mat) < 3000){
    enableFast = FALSE
    message("The number of cells in your dataset is less than 3,000. Fast mode has been disabled.")
    } else {
  message("The number of cells in your dataset exceeds 3,000. CytoTRACE will now be run in fast mode (see documentation). You can multi-thread this run using the 'ncores' flag. To disable fast mode, please indicate 'enableFast = FALSE'.")
  }
  #Checkpoint: NAs and poor quality genes
  pqgenes <- is.na(rowSums(mat>0)) | apply(mat, 1, var) == 0
  num_pqgenes <- length(which(pqgenes == TRUE))
  mat <- mat[!pqgenes,]
  if(num_pqgenes>0){
    warning(paste(num_pqgenes, "genes have zero expression in the matrix and were filtered"))
  }

  #Subsample routine
  if(enableFast == FALSE){
    size <- ncol(mat)
  } else if (enableFast == TRUE & subsamplesize < ncol(mat)){
    size <- subsamplesize
  } else if (enableFast == TRUE & subsamplesize >= ncol(mat)){
    stop("Please choose a subsample size less than the number of cells in dataset.")
  }

  chunk <- round(ncol(mat)/size)
  subsamples <- split(1:ncol(mat), sample(factor(1:ncol(mat) %% chunk)))
  message(paste("CytoTRACE will be run on", chunk, "sub-sample(s) of approximately",
                round(mean(unlist(lapply(subsamples, length)))), "cells each using", min(chunk, ncores),"/", ncores, "core(s)"))

  message(paste("Pre-processing data and generating similarity matrix..."))
  batches <- parallel::mclapply(subsamples, mc.cores = min(chunk, ncores), function(subsample){
    #Checkpoint: log2-normalization
    mat <- mat[,subsample]
    batch <- batch[subsample]

    if(max(mat)<50){
      mat <- 2^mat - 1
    }

    #Checkpoint: ERCC standards
    if(length(grep("ERCC-", rownames(mat)))>0){
      mat <- mat[-grep("ERCC-", rownames(mat)),]
    }

    #Checkpoint: Sequencing depth normalization
    mat <- t(t(mat)/apply(mat, 2, sum))*1000000

    #Checkpoint: NAs and poor quality cells
    pqcells <- is.na(apply(mat>0, 2, sum)) | apply(mat>0, 2, sum) <= 10
    num_pqcells <- length(which(pqcells == TRUE))
    mat <- mat[,!pqcells]

    #Checkpoint: log2-normalize
    mat <- log(mat+1,2)
    mat <- data.matrix(mat)

    #Calculate pre-batch corrected gene counts
    counts <- apply(mat>0, 2, sum)
    #Checkpoint: Batch correction
    if(ncol(a1) == length(a2)){
      #filter poor quality cells from batch vector
      batch <- batch[!pqcells]

      #Run Combat
      suppressMessages(mat <- sva::ComBat(mat, batch, c()))
      mat <- data.matrix(mat)

      #Replace negative values after batch correction with zeroes for compatibility with downstream steps
      mat[which(mat<0)] <- 0
    }
    #Rescale each single cell with gene counts to convert relative transcript abundances to absolute RNA content prior to cell lysis (credit: Census, Qiu et al., 2017)
    census_normalize <- function(mat, counts) {
      xnl <- 2^data.matrix(mat) - 1
      rs <- apply(xnl, 2, sum)
      rnorm <- t(t(xnl) * counts/rs)
      A <- log(rnorm+1,2)
      return(A)
    }

    mat2 <- census_normalize(mat, counts)
    #Function to identify the most variable genes
    mvg <- function(matn) {
      A <- matn
      n_expr <- rowSums(A > 0);
      A_filt <- A[n_expr >= 0.05 * ncol(A),];
      vars <- apply(A_filt, 1, var);
      means <- apply(A_filt, 1, mean);
      disp <- vars / means;
      last_disp <- tail(sort(disp), 1000)[1];
      A_filt <- A_filt[disp >= last_disp,];

      return(A_filt)
    }

    #Filter out cells not expressing any of the 1000 most variable genes
    mat2.mvg <- mvg(mat2)
    rm1 <- colSums(mat2.mvg) == 0
    mat2 <- mat2[, !rm1]
    counts <- counts[!rm1]

    #Calculate similarity matrix
    similarity_matrix_cleaned <- function(similarity_matrix){
      D <- similarity_matrix
      cutoff <- mean(as.vector(D))
      diag(D) <- 0;
      D[which(D < 0)] <- 0;
      D[which(D <= cutoff)] <- 0;
      Ds <- D
      D <- D / rowSums(D);
      D[which(rowSums(Ds)==0),] <- 0
      return(D)
    }
    D <- similarity_matrix_cleaned(HiClimR::fastCor(mvg(mat2)))

    return(list(mat2 = mat2,counts = counts, D = D))
  }
  )
  #Prepare for downstream steps
  mat2 <- do.call(cbind, lapply(batches, function(x) x$mat2))
  counts <- do.call(c, lapply(batches, function(x) x$counts))
  filter <- colnames(a1)[-which(colnames(a1) %in% colnames(mat2))]
  if(length(filter)>0){
    warning(paste(length(filter), "poor quality cells were filtered based on low or no expression. See 'filteredCells' in returned object for names of filtered cells."))
  }
  #Calculate gene counts signature (GCS) or the genes most correlated with gene counts
  message("Calculating gene counts signature...")
  ds2 <- sapply(1:nrow(mat2), function(x) ccaPP::corPearson(mat2[x,],counts))
  names(ds2) <- rownames(mat2)
  gcs <- apply(mat2[which(rownames(mat2) %in% names(rev(sort(ds2))[1:200])),],2,mean)

  samplesize <- unlist(lapply(lapply(batches, function(x) x$counts), length))
  gcs2 <- split(gcs, as.numeric(rep(names(samplesize), samplesize)))
  D2 <- lapply(batches, function(x) x$D)

  #Regress gene counts signature (GCS) onto similarity matrix
  regressed <- function(similarity_matrix_cleaned, score){
    out <- nnls::nnls(similarity_matrix_cleaned,score)
    score_regressed <- similarity_matrix_cleaned %*% out$x
    return(score_regressed)
  }

  #Apply diffusion to regressed GCS using similarity matrix
  diffused <- function(similarity_matrix_cleaned, score, ALPHA = 0.9){
    vals <- score
    v_prev <- rep(vals);
    v_curr <- rep(vals);

    for(i in 1:10000) {
      v_prev <- rep(v_curr);
      v_curr <- ALPHA * (similarity_matrix_cleaned %*% v_curr) + (1 - ALPHA) * vals;

      diff <- mean(abs(v_curr - v_prev));
      if(diff <= 1e-6) {
        break;
      }
    }
    return(v_curr)
  }

  message("Smoothing values with NNLS regression and diffusion...")
  cytotrace <- parallel::mclapply(1:length(D2), mc.cores = ncores, function(i) {
    gcs_regressed <- regressed(D2[[i]], gcs2[[i]])
    gcs_diffused <- diffused(D2[[i]], gcs_regressed)
    cytotrace <- rank(gcs_diffused)
  }
  )

  cytotrace <- cytotrace_ranked <- unlist(cytotrace)
  cytotrace <- range01(cytotrace)

  #Calculate genes associated with CytoTRACE
  cytogenes <- sapply(1:nrow(mat2),
                          function(x) ccaPP::corPearson(mat2[x,], cytotrace))
  names(cytogenes) <- rownames(mat2)
  message("Calculating genes associated with CytoTRACE...")

  #Final steps
  names(cytotrace) <- names(gcs) <- names(counts) <- colnames(mat2)
  cytotrace <- cytotrace[colnames(a1)]; cytotrace_ranked <- cytotrace_ranked[colnames(a1)]; gcs <- gcs[colnames(a1)]; counts <- counts[colnames(a1)]

  mat2 <- t(data.frame(t(mat2))[colnames(a1),])
  names(cytotrace) <- names(cytotrace_ranked) <- names(gcs) <- names(counts) <- colnames(mat2) <- colnames(a1)

  message("Done")
  return(list(CytoTRACE = cytotrace, CytoTRACErank = cytotrace_ranked, cytoGenes = sort(cytogenes, decreasing = T), GCS = gcs, gcsGenes = sort(ds2, decreasing = T),
              Counts = counts, filteredCells = filter, exprMatrix = mat2))
}


####get proportions
get_proportion=function(seurat=NULL,cluster_identity='seurat_clusters', sample_identity='condition',condition.1=NULL, condition.2=NULL,display=TRUE){
    library("scProportionTest")
    prop_test <- sc_utils(seurat)
    
    prop_test <- permutation_test(
    prop_test, cluster_identity = cluster_identity,
    sample_1 = condition.1, sample_2 = condition.2,
    sample_identity = sample_identity)
    
    p=as.data.frame(prop_test@results)
    
    if(display ==TRUE){
        pl=permutation_plot(prop_test)
        plot(pl)
    }
        
    seurat[['diff_proportions']]=seurat[[cluster_identity]][,cluster_identity]%in%c(subset(p, subset = permutation.FDR<0.05 & abs(permutation.obs_log2FD)>0.58 & permutation.pval<0.05)[,'permutation.clusters'])
    seurat[[paste0('more_in_',condition.1)]]=seurat[[cluster_identity]][,cluster_identity]%in%c(subset(p, subset = permutation.FDR<0.05 & permutation.obs_log2FD < -0.58)[,'permutation.clusters'])
    seurat[[paste0('more_in_',condition.2)]]=seurat[[cluster_identity]][,cluster_identity]%in%c(subset(p, subset = permutation.FDR<0.05 & permutation.obs_log2FD > 0.58)[,'permutation.clusters'])
        
        return(seurat)
        

    
}
                      
                      ####get proportions Scanpy
get_proportion_scanpy=function(seurat=NULL,cluster_identity='seurat_clusters', sample_identity='condition',condition.1=NULL, condition.2=NULL,display=TRUE){
    library("scProportionTest")
    prop_test <- sc_utils(seurat)
    
    prop_test <- permutation_test(
    prop_test, cluster_identity = cluster_identity,
    sample_1 = condition.1, sample_2 = condition.2,
    sample_identity = sample_identity)
    
    p=as.data.frame(prop_test@results)
    
    if(display ==TRUE){
        pl=permutation_plot(prop_test)
        plot(pl)
    }
        
    seurat[['diff_proportions']]=seurat[[cluster_identity]][,cluster_identity]%in%c(subset(p, subset = permutation.FDR<0.05 & abs(permutation.obs_log2FD)>0.58 & permutation.pval<0.05)[,'permutation.clusters'])
    seurat[[paste0('more_in_',condition.1)]]=seurat[[cluster_identity]][,cluster_identity]%in%c(subset(p, subset = permutation.FDR<0.05 & permutation.obs_log2FD < -0.58)[,'permutation.clusters'])
    seurat[[paste0('more_in_',condition.2)]]=seurat[[cluster_identity]][,cluster_identity]%in%c(subset(p, subset = permutation.FDR<0.05 & permutation.obs_log2FD > 0.58)[,'permutation.clusters'])
        
        return(seurat)
        

    
}
                      
                      
##Adapted from PAgoda2 to seurat

makeGeneKnnGraph = function(seurat=NULL,nPcs=100, center=TRUE, fastpath=TRUE, maxit=1000, k=30, n.cores=1, verbose=TRUE) {
    library(irlba)   
    
       x <- seurat@assays$RNA@counts
    
    

      # TODO: factor out gene PCA calculation
      # Do the PCA
      nPcs <- min(nrow(x)-1,ncol(x)-1,nPcs)
      if (center) {
          cm <- Matrix::colMeans(x)
          pcs <- irlba(x, nv=nPcs, nu =0, center=cm, right_only = FALSE, fastpath = fastpath, maxit= maxit, reorth = TRUE)
      } else {
         pcs <- irlba(x, nv=nPcs, nu =0, right_only = FALSE, fastpath = fastpath, maxit= maxit, reorth = TRUE)
      }
      rownames(pcs$v) <- colnames(x)

      # Optional centering
      if (center) {
        pcs$center <- cm
        pcas <- as.matrix(t(t(x %*% pcs$v) - t(cm  %*% pcs$v)))
      } else {
        pcas <- as.matrix(x %*% pcs$v)
      }

      # Keep the names
      rownames(pcas) <- rownames(x)
      colnames(pcas) <- paste0('PC',seq(ncol(pcas)))

      # Save into genegraphs slot
      #genegraphs$genePCs <- pcs
      #genegraphs$geneRotated <- pcas

      # Using cosine distance only here
      if (center) {
        pcas <- pcas - Matrix::rowMeans(pcas)
      }
      xn <- N2R::Knn(pcas, k, nThreads= n.cores, verbose=verbose)
      diag(xn) <- 0 # Remove self edges
      xn <- as(xn,'dgTMatrix') # will drop 0s
      # Turn into a dataframe, convert from correlation distance into weight
      df <- data.frame('from'=rownames(pcas)[xn@i+1],'to'=rownames(pcas)[xn@j+1],'w'=pmax(1-xn@x,0),stringsAsFactors=FALSE)

      
    return(df)
    }



#####Proportion differences between 2 conditions

get_proportion=function(seurat=NULL,cluster_identity='seurat_clusters', sample_identity='condition',condition.1=NULL, condition.2=NULL,display=TRUE){
    library("scProportionTest")
    prop_test <- sc_utils(seurat)
    
    prop_test <- permutation_test(
    prop_test, cluster_identity = cluster_identity,
    sample_1 = condition.1, sample_2 = condition.2,
    sample_identity = sample_identity)
    
    p=as.data.frame(prop_test@results)
    
    if(display ==TRUE){
        pl=permutation_plot(prop_test)
        plot(pl)
    }
        
    seurat[['diff_proportions']]=seurat[[cluster_identity]][,cluster_identity]%in%c(subset(p, subset = permutation.FDR<0.05 & abs(permutation.obs_log2FD)>0.58 & permutation.pval<0.05)[,'permutation.clusters'])
    seurat[[paste0('more_in_',condition.1)]]=seurat[[cluster_identity]][,cluster_identity]%in%c(subset(p, subset = permutation.FDR<0.05 & permutation.obs_log2FD < -0.58)[,'permutation.clusters'])
    seurat[[paste0('more_in_',condition.2)]]=seurat[[cluster_identity]][,cluster_identity]%in%c(subset(p, subset = permutation.FDR<0.05 & permutation.obs_log2FD > 0.58)[,'permutation.clusters'])
        
        return(seurat)
        

    
}