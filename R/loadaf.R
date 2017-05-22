# --------------------------------------
# The Sync-class and related functions -
#---------------------------------------

# class definition -----
setClass(Class="sync",

         # member variables
         slots=c(
           gen="numeric",
           repl="numeric",
           isAF="logical",
           alleles="data.table"),

         # default constructor
         prototype=list(
           gen=numeric(0),
           repl=numeric(0),
           isAF=logical(0),
           alleles=NULL),

         # validation
         validity=function(object) {
           # length of 'gen', 'repl' and 'isAF' do not all match
           if(!(length(object@gen) == length(object@repl) && length(object@gen) == length(object@isAF))) {
             return(paste0("Lengths of 'gen' (", length(object@gen), "), 'repl' (", length(object@repl), ") and 'isAF' (", length(object@isAF), ") do not all match."))
           }
           # first six column names are not set as expected
           if(!is.null(object@alleles) && !isTRUE(all.equal(colnames(object@alleles)[1:6], c("chr", "pos", "ref", "major", "minor", "rising")))) {
             return(paste0("The names of the first five columns in 'alleles' must be 'c(\"chr\", \"pos\", \"ref\", \"major\", \"minor\", \"rising\") in exactly that order"))
           }
           # column number of 'alleles' does not match the length of 'gen'
           if((is.null(object@alleles) && length(object@gen) != 0) || (!is.null(object@alleles) && ncol(object@alleles) != length(object@gen)+6)) {
             return(paste0("Column number of 'alleles' (", ncol(object@alleles), ") has to be equal to the length of 'gen' (", length(object@gen), ") + 6."))
           }
           # neither allele frequencies nor sequence coverages are provided
           if(length(object@gen) == 0) {
             return("Neither allele frequencies nor sequence coverage are provided")
           }
           # data types of columns in 'alleles' do not match
           for(cc in c(2, seq(7, length(object@gen)+6))) {
             if(!is.numeric(object@alleles[[cc]])) {
               return(paste0("Column ", cc, " in 'alleles' should be numeric."))
             }
           }
           for(cc in c(1, 3:6)) {
             if(!is.character(object@alleles[[cc]])) {
               return(paste0("Column ", cc, " in 'alleles' should be character."))
             }
           }

           # allele frequency is outside the range [0, 1]
           for(cc in which(object@isAF)+6) {
             if(!all(object@alleles[[cc]] >= 0, object@alleles[[cc]] <= 1, na.rm=TRUE)) {
               return(paste0("Allele frequencies in column ", cc, " are outside the range [0, 1]."))
             }
           }

           # coverage is < 0
           for(cc in which(!object@isAF)+6) {
             if(!all(object@alleles[[cc]] >= 0, na.rm=TRUE)) {
               return(paste0("Sequence coverages in column ", cc, " are < 0."))
             }
           }

           # object is valid
           return(TRUE)
         }
)

# initializer -----
setMethod("initialize",
          signature="sync",
          definition=function(.Object, gen, repl, isAF, alleles) {
            .Object@gen <- gen
            .Object@repl <- repl
            .Object@isAF <- isAF
            .Object@alleles <- alleles
            # call the inspector
            validObject(.Object)
            # if 'alleles' is specified then add column with 'posID' and set key
            if(!is.null(.Object@alleles)) {
              .Object@alleles[,posID:=paste(chr, pos, sep=".")]
              setkey(.Object@alleles, posID)
            }
            return(.Object)
          }
)


# merge sync-objects ----
setMethod("merge",
          signature=c(x="sync", y="sync"),
          definition=function(x, y, ...) {

          } )

# related methods -----
is.sync <- function(x) {
  return(inherits(x, "sync"))
}



af.traj <- function(sync, chr, pos, repl) {
  # if 'sync' is not inherited from class 'sync' then stop execution
  if(!is.sync(sync))
    stop("Argument 'sync' is not a sync-object.")

  # extract selected SNPs
  sel <- paste(chr <- as.character(chr), pos <- as.numeric(pos), sep=".")
  subAlleles <- sync@alleles[sel]

  # extract time series data for specific replicate
  traj.mat <- function(r, include=FALSE) {
    cols <- which(sync@repl == r & sync@isAF)
    res <- subAlleles[,cols[order(sync@gen[cols])]+6,with=FALSE]
    res <- as.matrix(setnames(res, sub("\\.freq", "", colnames(res))))
    rownames(res) <- if(include) paste0(sel, ".R", r) else sel
    colnames(res) <- sub("\\.R[0-9]+", "", colnames(res))
    return(res)
  }

  # if only one replicate is specified then create 1xt matrix for trajectory with t time points
  if(length(repl <- as.numeric(repl)) == 1) {
    traj <- traj.mat(repl, include=FALSE)
    # if multiple replicates are specified
  } else if(length(repl) > 1) {
    # then create list of trajectory matrices, one for each replicate
    traj <- foreach(r=repl) %do% {
      traj.mat(r, include=TRUE)
    }
    # if the number of time points (columns) is equal for all replicates then combine all matrices to one
    if(all((x <- sapply(traj, ncol)) == x[1])) {
      traj <- Reduce(rbind, traj)
    }
  }

  return(traj)
}

af <- function(sync, repl, gen) {
  # if 'sync' is not inherited from class 'sync' then stop execution
  if(!is.sync(sync))
    stop("Argument 'sync' is not a sync-object.")

  # extract allele frequencies of specified replicates and time points
  cols <- which(sync@gen %in% (gen <- as.numeric(gen)) & sync@repl %in% (repl <- as.numeric(repl)) & sync@isAF)
  if(length(cols) == 0)
    stop("The combination of 'repl' (", paste(repl, sep=","), ") and 'gen' (", paste(gen), ") does not match the data available in 'sync'.")
  subAlleles <- sync@alleles[,cols[order(sync@repl[cols], sync@gen[cols])]+6,with=FALSE]
  afMat <- as.matrix(subAlleles)
  rownames(afMat) <- sync@alleles$posID

  return(afMat[order(sync@alleles$chr, sync@alleles$pos),])
}

coverage <- function(sync, repl, gen) {
  # if 'sync' is not inherited from class 'sync' then stop execution
  if(!is.sync(sync))
    stop("Argument 'sync' is not a sync-object.")

  # extract allele frequencies of specified replicates and time points
  cols <- which(sync@gen %in% (gen <- as.numeric(gen)) & sync@repl %in% (repl <- as.numeric(repl)) & !sync@isAF)
  if(length(cols) == 0)
    stop("The combination of 'repl' (", paste(repl, sep=","), ") and 'gen' (", paste(gen), ") does not match the data available in 'sync'.")
  subAlleles <- sync@alleles[,cols[order(sync@repl[cols], sync@gen[cols])]+6,with=FALSE]
  afMat <- as.matrix(subAlleles)
  rownames(afMat) <- sync@alleles$posID

  return(afMat[order(sync@alleles$chr, sync@alleles$pos),])
}

alleles <- function(sync) {
  # if 'sync' is not inherited from class 'sync' then stop execution
  if(!is.sync(sync))
    stop("Argument 'sync' is not a sync-object.")

  return(sync@alleles[order(chr, pos),1:6,with=FALSE])
}

splitLocusID <- function(id, sep=".") {
  splitMat <- stri_split(id <- as.character(id), fixed=(sep <- as.character(sep)), simplify=TRUE)

  if(!is.matrix(splitMat) || dim(splitMat)[1] != length(id) || dim(splitMat)[2] != 2)
    stop("Locus IDs could not be split successfully.")

  return(data.table(chr=splitMat[,1], pos=as.numeric(splitMat[,2])))
}


# ---------------------------------------------------
# Convert Sync-file to allele frequency data object -
# ---------------------------------------------------

read.sync <- function(file, gen, repl, rising=FALSE) {

  cat("Reading sync file ...\n")
  # load sync-file
  syncDt <- read.delim(file, sep="\t", header=FALSE, stringsAsFactors=FALSE)
  setDT(syncDt)
  # if either 'gen' or 'repl' are not of propper length then stop
  if((ncol(syncDt)-3) %% length(gen) != 0 || (ncol(syncDt)-3) %% length(repl) != 0)
    stop("Either 'gen' (", length(gen), ") or 'repl' (", length(repl), ") is not a multiple of the number of populations (", ncol(syncDt)-3, ") specified in the sync-file.")
  gc()

  cat("Extracting biallelic counts ...\n")
  # extract numeric allele counts for A:T:C:G
  syncCnts <- lapply(syncDt[,-1:-3,with=FALSE], function(col) {
    matrix(as.numeric(stri_split(col, fixed=":", simplify=TRUE)), ncol=6)[,1:4]
  } ) # <- bottleneck

  # get rid of data that is no longer needed
  snpCnt <- nrow(syncDt)
  chr <- as.character(syncDt$V1)
  pos <- syncDt$V2
  ref <- syncDt$V3
  popCnt <- ncol(syncDt)-3
  rm(syncDt)
  gc()

  # sum up counts across populations (time points and replicates) and add e (random uniform >=0 & <= 0.99) to make each count value unique
  sumCnts <- Reduce("+", syncCnts)
  sumCnts <- sumCnts + runif(nrow(sumCnts)*ncol(sumCnts), min=0, max=0.99)
  # deterime allele ranks for each position
  alleleRank <- rowRanks(sumCnts, ties.method="min")
  # extract 2 most common alleles (considering all populations)
  alleleCnts <- lapply(syncCnts, function(pop) {
    cbind(major=t(pop)[t(alleleRank) == 4], minor=t(pop)[t(alleleRank) == 3])
  } ) # <- bottleneck

  cat("Creating result object ...\n")
  # compute chromosome IDs (to later replace character vector by numeric one)
  chrNames <- unique(chr)
  chrID <- 1:length(chrNames)
  names(chrID) <- chrNames

  # extract major and minor allele for each position
  syncCntCol <- 1:4
  names(syncCntCol) <- c("A", "T", "C", "G")
  alleles <- data.table(chr=chr, pos=pos, ref=ref,
                        major=names(syncCntCol)[which(t(alleleRank) == 4) - 4*seq(0, nrow(alleleRank)-1)],
                        minor=names(syncCntCol)[which(t(alleleRank) == 3) - 4*seq(0, nrow(alleleRank)-1)],
                        rising=NA_character_)

  # combine generation and replicate info
  popInfo <- data.table(pop=1:popCnt, gen=gen, repl=repl)

  # add allele frequency and sequence coverage for each population
  for(r in unique(repl)) {
    for(i in seq(1:nrow(popInfo))[popInfo$repl == r]) {
      seqCov <- rowSums(alleleCnts[[i]])
      alleles[,paste0("F", popInfo$gen[i], ".R", r, ".freq"):=alleleCnts[[i]][,"minor"]/seqCov]
      alleles[,paste0("F", popInfo$gen[i], ".R", r, ".cov"):=seqCov]
    }
  }

  # if required then polarize allele counts for the rising allele
  if(rising && length(unique(popInfo$gen)) > 1) {
    ugens <- unique(popInfo$gen)
    minGen <- min(popInfo$gen)

    # if minGen allele frequency column is not available for all replicates then stop execution
    if(sum(grepl(paste0("F", minGen, "\\.[R0-9]+\\.freq"), colnames(alleles))) != length(unique(repl)))
      stop("Not all replicates provide allele frequency estimates at generation F", minGen)

    # calculate mean allele frequency change per SNP and replicate
    meanAF <- foreach(r=unique(repl), .combine=cbind, .final=function(x) { if(is.matrix(x)) return(rowMeans(x)) else return(x) }) %do% {
      rowMeans(alleles[,grepl(paste0("F[0-9]+\\.R", r, "\\.freq"), colnames(alleles)),with=FALSE]-alleles[[which(grepl(paste0("F", minGen, "\\.[R0-9]+\\.freq"), colnames(alleles)))[1]]])
    }

    # polarize allele frequencies
    needsPolarization <- meanAF < 0
    for(pop in grep("F[0-9]+\\.R[0-9]+\\.freq", colnames(alleles), value=TRUE)) {
      alleles[,eval(pop):=ifelse(needsPolarization, 1-alleles[[pop]], alleles[[pop]])]
    }

    # set column with rising allele
    alleles[,rising:=ifelse(needsPolarization, alleles$major, alleles$minor)]
  }

  # return sync-object for loaded sync-file
  return(new(Class="sync",
             gen=as.numeric(sub("F([0-9]+)\\.R[0-9]+.*", "\\1", colnames(alleles)[-1:-6])),
             repl=as.numeric(sub(".*\\.R([0-9]+)\\..*", "\\1", colnames(alleles)[-1:-6])),
             isAF=grepl(".*\\.freq$", colnames(alleles)[-1:-6]),
             alleles=alleles))
}
