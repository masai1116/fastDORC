## fastDORC funtion for fast identifying peak-gene association from SHARE-seq or other scATAC-RNA datasets
## for related question, email vinay_kartha@g.harvard.edu or sai@broadinstitute.org

library(SummarizedExperiment)
library(Matrix)
hg19TSSRanges <- readRDS("hg19TSSRanges.rds")
hg38TSSRanges <- readRDS("hg38TSSRanges.rds")
mm10TSSRanges <- readRDS("mm10TSSRanges.rds")
hg19.ncRNATSSRanges <- readRDS("hg19.ncRNATSSRanges.rds")
atac.se.exmaple <- readRDS("atac.se.example.rds")
RNA.example <- readRDS("RNA.example.rds")
# RNA and atac objects should have same number of cells

fastGenePeakcorr <- function(ATAC.se,
                             RNAmat, 
                             genome, # Must be one of "hg19", "mm10", "hg38", or "hg19.ncRNA"; hg19.ncRNA contains both mRNA and ncRNA annotation
                             geneList=NULL, # 2 or more valid gene symbols
                             windowPadSize=50000,
                             normalizeATACmat=TRUE,
                             nCores=6,
                             keepMultiMappingPeaks=FALSE,
                             n_bg=100,
                             p.cut=NULL # Optional, if specified, will only return sig hits
) {
  stopifnot(inherits(ATAC.se,"RangedSummarizedExperiment"))
  stopifnot(inherits(RNAmat,c("Matrix","matrix")))
  if(!all.equal(ncol(ATAC.se),ncol(RNAmat)))
    stop("Input ATAC and RNA objects must have same number of cells")
  # Function needs rownames for both matrices or gives error
  rownames(ATAC.se) <- paste0("Peak",1:nrow(ATAC.se))
  peakRanges.OG <- granges(ATAC.se)
  if(is.null(rownames(RNAmat)))
    stop("RNA matrix must have gene names as rownames")
  # Check for peaks/genes with 0 accessibility/expression
  if(any(Matrix::rowSums(assay(ATAC.se))==0)){
    message("Peaks with 0 accessibility across cells exist ..")
    message("Removing these peaks prior to running correlations ..")
    peaksToKeep <- Matrix::rowSums(assay(ATAC.se))!=0
    ATAC.se <- ATAC.se[peaksToKeep,] # Subset ranges
    message("Important: peak indices in returned gene-peak maps are relative to original input SE")
  }
  ATACmat <- assay(ATAC.se) # Rownames preserved
  # Normalize peak counts
  if (normalizeATACmat) 
    ATACmat <- centerCounts(ATACmat)
    peakRanges <- granges(ATAC.se) # Peak ranges
  if(any(Matrix::rowSums(RNAmat)==0)){
    message("Genes with 0 expression across cells exist ..")
    message("Removing these genes prior to running correlations ..")
    genesToKeep <- Matrix::rowSums(RNAmat)!=0
    RNAmat <- RNAmat[genesToKeep,]
  }
  cat("Number of peaks in ATAC data:",nrow(ATACmat),"\n")
  cat("Number of genes in RNA data:",nrow(RNAmat),"\n")
  if (!genome %in% c("hg19", "hg38", "mm10", "hg19.ncRNA")) 
    stop("You must specify one of hg19, hg38, hg19.ncRNA or mm10 as a genome build for currently supported TSS annotations..\n")
  switch(genome, hg19 = {
    TSSg <- hg19TSSRanges
  }, hg38 = {
    TSSg <- hg38TSSRanges
  }, mm10 = {
    TSSg <- mm10TSSRanges
  }, hg19.ncRNA = {
    TSSg <- hg19.ncRNATSSRanges
  })
  # Keep genes that have annotation and are in RNA matrix
  names(TSSg) <- as.character(TSSg$gene_name)
  if(!is.null(geneList)){
    if(length(geneList)==1)
      stop("Please specify more than 1 valid gene symbol")
    if(any(!geneList %in% names(TSSg))){
      cat("One or more of the gene names supplied is not present in the TSS annotation specified: \n")
      cat(geneList[!geneList %in% names(TSSg)], sep = ", ")
      cat("\n")
      stop()
    }
    TSSg <- TSSg[geneList]
  }
  # Checking in case some genes in RNA don't overlap our TSS annotations
  genesToKeep <- intersect(names(TSSg),rownames(RNAmat))
  cat("\nNum genes overlapping TSS annotation and RNA matrix being considered: ",length(genesToKeep),"\n")
  # Match gene order in RNA matrix and TSS ranges
  RNAmat <- RNAmat[genesToKeep,]
  TSSg <- TSSg[genesToKeep]
  # Pad TSS by this much *either side*
  TSSflank <- GenomicRanges::flank(TSSg, 
                                   width = windowPadSize,
                                   both = TRUE)
  # Get peak summit
  cat("\nTaking peak summits from peak windows ..\n")
  peakSummits <- resize(peakRanges,width = 1,fix = "center")
  # Find overlap of all peaks to all genes given window
  # Subject is Peaks, query is Gene
  cat("Finding overlapping peak-gene pairs ..\n")
  genePeakOv <- findOverlaps(query = TSSflank,subject = peakSummits)
  numPairs <- length(genePeakOv)
  cat("Found ",numPairs,"total gene-peak pairs for given TSS window ..\n")
  cat("Number of peak summits that overlap any gene TSS window: ",length(unique(subjectHits(genePeakOv))),"\n")
  cat("Number of gene TSS windows that overlap any peak summit: ",length(unique(queryHits(genePeakOv))),"\n\n")
  # For each gene, determine observed correlation of each overlapping peak to its associated gene (gene expression)
  # For each of those genes, also determine correlation based on background peaks (run in parallel) and save
  # Do this in parallel, and get data frame of gene-peak-pearson values
  # Fetch background peaks for each peak tested (i.e. that has overlap in window with gene)
  cat("Determining background gene-peak pairs ..\n")
  if(is.null(rowData(ATAC.se)$bias)){
    # Get GC content if not pre-determined in SE
    if(genome %in% "hg19")
      myGenome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
    if(genome %in% "mm10")
      myGenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
    if(genome %in% "hg38")
      myGenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    if(genome %in% "hg19.ncRNA")
      myGenome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
    ATAC.se <- chromVAR::addGCBias(ATAC.se,genome=myGenome) 
  }
  cat("Using ",n_bg," iterations ..\n\n")
  start.time <- Sys.time()
  # Randomly select background peak-gene pairs
  # Number of (reference) gene-peak pairs to precompute correlation scaffolds for
  numBgPairs <- 100000
  #numBgPairs <- 300000
  set.seed(123)
  bgPairGenes <- sample(1:length(genesToKeep), numBgPairs, replace = TRUE)
  bgPairPeaks <- sample(1:length(peakSummits), numBgPairs, replace = TRUE)
  # Get GC content, accessibilities of background peaks. Also get expression level of background genes.
  bgPairFeatures <- data.frame(GC = ATAC.se@rowRanges$bias[bgPairPeaks],
                               accessibility = rowMeans(ATACmat)[bgPairPeaks],
                               expression = rowMeans(RNAmat)[bgPairGenes])
  # Also get GC content, accessibilities, expression levels for observed peak-gene pairs
  obPairGenes <- genePeakOv@from
  obPairPeaks <- genePeakOv@to
  obPairFeatures <- data.frame(GC = ATAC.se@rowRanges$bias[obPairPeaks],
                               accessibility = rowMeans(ATACmat)[obPairPeaks],
                               expression = rowMeans(RNAmat)[obPairGenes])
  # Rescale these features so GC/accessibility/expression roughly have the same weight
  allPairFeatures <- scale(rbind(as.matrix(bgPairFeatures), as.matrix(obPairFeatures)))
  bgPairFeatures <- allPairFeatures[1:nrow(bgPairFeatures),]
  obPairFeatures <- allPairFeatures[(nrow(bgPairFeatures) + 1):
                                      (nrow(bgPairFeatures) + nrow(obPairFeatures)),]
  # Find all background pairs of the observed pairs by searching for nearest neighbor pairs in GC-accessibility-expression space
  bgPairsInds <- FNN::get.knnx(data = bgPairFeatures, query = obPairFeatures, k=n_bg)$nn.index
  metric <- "spearman"
  pairsPerChunk <- 500
  # Compute the correlation for all background pairs
  pairCorrs <- list()
  for(pairs in c("bg", "ob")){
    # Divide the list of peak-gene pairs in to chunks
    if(pairs == "bg") {
      numPairs <- length(bgPairGenes)
      cat("\nComputing background peak-gene pair correlations ..\n")
    } else if(pairs == "ob") {
      numPairs <- length(obPairGenes)
      cat("\nComputing observed peak-gene pair correlations ..\n")
    }
    starts <- seq(1, numPairs, pairsPerChunk)
    ends <- starts  + pairsPerChunk -1
    ends[length(ends)] <- numPairs
    chunkList <- mapply(c, starts, ends, SIMPLIFY = FALSE)
    # Parallelized calculation of peak-gene correlation for each chunk
    if(pairs == "bg") {
      corPairs <- data.frame(Gene = bgPairGenes, Peak = bgPairPeaks,stringsAsFactors = FALSE)
    } else if(pairs == "ob") {
      corPairs <- data.frame(Gene = obPairGenes, Peak = obPairPeaks,stringsAsFactors = FALSE)
    }
    library(doParallel)
    if(nCores > 1)
      cat("Running in parallel using ",nCores," cores ..\n\n")
    opts <- list()
    pb <- txtProgressBar(min = 0, max = length(chunkList), style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    time_elapsed <- Sys.time()
    cl <- parallel::makeCluster(nCores)
    clusterEvalQ(cl, .libPaths())
    doSNOW::registerDoSNOW(cl)
    corList <- foreach(x=1:length(chunkList),
                       .options.snow = opts, 
                       .export = c(".chunkCore"),
                       .packages = c("dplyr","Matrix")) %dopar%  {
                         corrs <- .chunkCore(chunk=chunkList[[x]],A=ATACmat,R=RNAmat,O=corPairs,met=metric)
                         return(corrs)
                       }
    stopCluster(cl)
    if(any(unlist(sapply(corList,is.null)))){
      message("One or more of the chunk processes failed unexpectedly (returned NULL) ..")
      message("Please check to see you have enough cores/memory allocated")
      message("Also make sure you have filtered down to non-zero peaks/genes")
    }
    pairCorrs[[pairs]] <- unlist(corList)
  }
  # Get correlation pvalues by comparing observes correlation to background distribution
  cat("\nCalculating correlation P-values based on background distribution ..\n")
  pvals <- pbmcapply::pbmcmapply(
    function(pair_ind){
      obcor <- pairCorrs[["ob"]][pair_ind]
      bgCorrs <- pairCorrs[["bg"]][bgPairsInds[pair_ind,]]
      pval <- 1-stats::pnorm(q = obcor, mean = mean(bgCorrs), sd = sd(bgCorrs))
    },
    1:length(obPairGenes),
    mc.cores = nCores
  )
  time_elapsed <- Sys.time() - start.time
  cat(paste("\nTime Elapsed: ", round(time_elapsed,2), units(time_elapsed), "\n"))
  corrResults <- data.frame(Peak = obPairPeaks,
                            Gene = obPairGenes,
                            rObs = pairCorrs[["ob"]],
                            pvalZ = pvals,
                            stringsAsFactors = FALSE)
  # Filter to positive correlations
  cat("Only considering positive associations ..\n")
  corrResults <- corrResults %>% filter(rObs > 0)
  if(!keepMultiMappingPeaks){
    # Remove multi-mapping peaks (force 1-1 mapping)
    cat("Keeping max correlation for multi-mapping peaks ..\n")
    corrResults <- corrResults %>% group_by(Peak) %>% filter(rObs==max(rObs))
  }
  # Swap gene number for gene symbol from TSS annotation lookup
  corrResults$Gene <- as.character(TSSg$gene_name)[corrResults$Gene]
  # Swap peak numbers to match reference input peak numbers
  # This only changes if some peaks had zero accessibility and were filtered out
  # Use rownames from reference
  corrResults$Peak <- as.numeric(BuenRTools::splitAndFetch(rownames(ATACmat)[corrResults$Peak],"Peak",2))
  cat("\nFinished!\n")
  # If there is a significance cutoff, only keep pairs that are significant
  if(!is.null(p.cut)){
    cat("Using significance cut-off of ",p.cut," to subset to resulting associations\n")
    corrResults <- corrResults[corrResults$pvalZ <= p.cut,] # Subset to significant correlations only
  }
  # Add ranges for peaks
  corrResults$PeakRanges <- paste(as.character(seqnames(peakRanges.OG[corrResults$Peak])),paste(start(peakRanges.OG[corrResults$Peak]),end(peakRanges.OG[corrResults$Peak]),sep="-"),sep=":")  
  return(corrResults %>% as.data.frame(stringsAsFactors=FALSE))
}

.chunkCore <- function(chunk,
                       A, # ATAC matrix
                       R, # RNA matrix
                       O, # Gene-Peak overlap pairing data.frame
                       met # Correlation method ("spearman" or "pearson")
){
  # Get indices of genes and peaks from overlap object for chunk
  # Assumes query hits are genes and subject hits are peaks in the overlap object
  geneIndices <- O$Gene[chunk[1]:chunk[2]]
  peakIndices <- O$Peak[chunk[1]:chunk[2]]
  pairnames <- cbind(rownames(A)[peakIndices],rownames(R)[geneIndices])
  uniquegenes <- unique(geneIndices)
  uniquepeaks <- unique(peakIndices)
  M1 <- as.matrix(t(A[uniquepeaks,,drop=FALSE])) # In case only 1 match, keep matrix structure
  M2 <- as.matrix(t(R[uniquegenes,,drop=FALSE])) # In case only 1 match, keep matrix structure
  # Peak x Gene correlation matrix, subset by peak-gene pair names to get corresponding correlation vector
  cor(x = M1,y = M2,method = met)[pairnames]
}

centerCounts <- function (obj, doInChunks = TRUE, chunkSize = 1000) 
{
  if (!class(obj) %in% c("SummarizedExperiment", "RangedSummarizedExperiment", 
                         "dgCMatrix", "dgeMatrix", "Matrix")) 
    stop("Supplied object must be either of class SummarizedExperiment or sparse Matrix ..\n")
  if (ncol(obj) > 10000) 
    doInChunks <- TRUE
  if (doInChunks) {
    cat("Centering counts for cells sequentially in groups of size ", 
        chunkSize, " ..\n\n")
    starts <- seq(1, ncol(obj), chunkSize)
  }
  else {
    starts <- 1
  }
  counts.l <- list()
  for (i in 1:length(starts)) {
    beginning <- starts[i]
    if (i == length(starts)) {
      ending <- ncol(obj)
    }
    else {
      ending <- starts[i] + chunkSize - 1
    }
    cat("Computing centered counts for cells: ", beginning, 
        " to ", ending, "..\n")
    if (class(obj) == "RangedSummarizedExperiment" | class(obj) == 
        "SummarizedExperiment") {
      m <- SummarizedExperiment::assay(obj[, beginning:ending])
    }
    else {
      m <- obj[, beginning:ending]
    }
    cellMeans <- Matrix::colMeans(m)
    cat("Computing centered counts per cell using mean reads in features ..\n\n")
    cCounts <- Matrix::t(Matrix::t(m)/cellMeans)
    counts.l[[i]] <- cCounts
    gc()
  }
  cat("Merging results..\n")
  centered.counts <- do.call("cbind", counts.l)
  cat("Done!\n")
  if (class(obj) == "RangedSummarizedExperiment" | class(obj) == 
      "SummarizedExperiment") {
    SummarizedExperiment::assay(obj) <- centered.counts
    return(obj)
  }
  else {
    return(centered.counts)
  }
}

cisCor <- fastGenePeakcorr(ATAC.se = atac.se.exmaple, 
                           RNAmat = RNA.example, 
                           genome = "hg19.ncRNA", 
                           # geneList = c("MIR4420", "LEF1", "GATA1"), 
                           windowPadSize = 500000, # recommend 50000 or 500000 on each side of TSSs
                           normalizeATACmat = TRUE,
                           nCores = 4,
                           keepMultiMappingPeaks = FALSE,
                           n_bg = 100,
                           p.cut = NULL
)
cisCor$PeakCoor <- atac.se.exmaple@rowRanges$peak[cisCor$Peak]
cisCor.sig <- cisCor[cisCor$pvalZ < 0.05, ]

# saveRDS(cisCor, "cisCor.rds")
# saveRDS(cisCor.sig, "cisCor.sig.rds")
