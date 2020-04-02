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
           polarization="character",
           alleles="data.table"),

         # default constructor
         prototype=list(
           gen=numeric(0),
           repl=numeric(0),
           isAF=logical(0),
           polarization=character(0),
           alleles=NULL),

         # validation
         validity=function(object) {
           # length of 'gen', 'repl' and 'isAF' do not all match
           if(!(length(object@gen) == length(object@repl) && length(object@gen) == length(object@isAF))) {
             return(paste0("Lengths of 'gen' (", length(object@gen), "), 'repl' (", length(object@repl), ") and 'isAF' (", length(object@isAF), ") do not all match."))
           }
           # check polarization
           if(!object@polarization %in% c("minor", "rising", "reference")) {
             return(paste0("Polarization needs to be either 'major', 'rising' or 'reference': ", object@polarization))
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
          definition=function(.Object, gen, repl, isAF, polarization, alleles) {
            .Object@gen <- gen
            .Object@repl <- repl
            .Object@isAF <- isAF
            .Object@polarization <- polarization
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

polarization <- function(sync) {
  sync@polarization
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

af <- function(sync, chr, pos, repl, gen) {
  # if 'sync' is not inherited from class 'sync' then stop execution
  if(!is.sync(sync))
    stop("Argument 'sync' is not a sync-object.")

  # determine for which replicate and time point allele frequencies should be extracted
  cols <- which(sync@gen %in% (gen <- as.numeric(gen)) & sync@repl %in% (repl <- as.numeric(repl)) & sync@isAF)
  if(length(cols) == 0)
    stop("The combination of 'repl' (", paste(repl, sep=","), ") and 'gen' (", paste(gen), ") does not match the data available in 'sync'.")

  # if required get AF for subset of all SNPs
  subAlleles <- sync@alleles
  snps <- NULL
  if(!missing(chr) && !missing(pos)) {
    snps <- paste(chr, pos, sep=".")
    subAlleles <- subAlleles[snps]
  }
  subAlleles <- subAlleles[,cols[order(sync@repl[cols], sync@gen[cols])]+6,with=FALSE]

  # convert data.table to matrix
  afMat <- as.matrix(subAlleles)
  rownames(afMat) <- if(is.null(snps)) sync@alleles$posID else snps

  # return allele frequencies ordered by chromosome and position, or as specified by 'chr' and 'pos' parameters
  return(if(is.null(snps)) afMat[order(sync@alleles$chr, sync@alleles$pos),] else afMat)
}

coverage <- function(sync, chr, pos, repl, gen) {
  # if 'sync' is not inherited from class 'sync' then stop execution
  if(!is.sync(sync))
    stop("Argument 'sync' is not a sync-object.")

  # determine for which replicate and time point sequence coverages should be extracted
  cols <- which(sync@gen %in% (gen <- as.numeric(gen)) & sync@repl %in% (repl <- as.numeric(repl)) & !sync@isAF)
  if(length(cols) == 0)
    stop("The combination of 'repl' (", paste(repl, sep=","), ") and 'gen' (", paste(gen), ") does not match the data available in 'sync'.")

  # if required get coverages for subset of all SNPs
  subAlleles <- sync@alleles
  snps <- NULL
  if(!missing(chr) && !missing(pos)) {
    snps <- paste(chr, pos, sep=".")
    subAlleles <- subAlleles[snps]
  }
  subAlleles <- subAlleles[,cols[order(sync@repl[cols], sync@gen[cols])]+6,with=FALSE]

  # convert data.table to matrix
  afMat <- as.matrix(subAlleles)
  rownames(afMat) <- if(is.null(snps)) sync@alleles$posID else snps

  # return coverages ordered by chromosome and position, or as specified by 'chr' and 'pos' parameters
  return(if(is.null(snps)) afMat[order(sync@alleles$chr, sync@alleles$pos),] else afMat)
}

alleles <- function(sync, chr, pos) {
  # if 'sync' is not inherited from class 'sync' then stop execution
  if(!is.sync(sync))
    stop("Argument 'sync' is not a sync-object.")

  # if 'chr' and 'pos' are specified then return only information for those loci
  if(!missing(chr) && !missing(pos)) {
    snps <- paste(chr, pos, sep=".")
    return(sync@alleles[snps][,1:6,with=FALSE])
  }

  # otherwise return information on all loci
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

read.sync <- function(file, gen, repl, polarization=c("minor", "rising", "reference"), keepOnlyBiallelic=FALSE) {

  polarization <- match.arg(polarization)

  cat("Reading sync file ...\n")
  # load sync-file
  syncDt <- fread(file, sep="\t", header=FALSE, stringsAsFactors=FALSE)
  setDT(syncDt)
  # if either 'gen' or 'repl' are not of propper length then stop
  if((ncol(syncDt)-3) %% length(gen) != 0 || (ncol(syncDt)-3) %% length(repl) != 0)
    stop("Either 'gen' (", length(gen), ") or 'repl' (", length(repl), ") is not a multiple of the number of populations (", ncol(syncDt)-3, ") specified in the sync-file.")

  cat("Extracting biallelic counts ...\n")

  cppFunction(plugins=c("cpp11"),{"NumericMatrix Sync2Cnts(CharacterVector sync) {
    NumericMatrix resMat(sync.length(), 4);

    std::string current;
    std::string token;
    size_t posBeg = 0, posEnd = 0;
    int cnt = 0;

    for(int i=0; i<sync.length(); i++) {

    cnt = 0;
    posBeg = 0;
    posEnd = 0;
    current = Rcpp::as<std::string>(sync(i));
    while ((posEnd = current.find(\":\", posBeg)) != std::string::npos && cnt <= 3) {
    token = current.substr(posBeg, posEnd);
    resMat(i, cnt) = std::stoi(token);
    posBeg = posEnd+1;
    cnt++;
    }
  }

    return resMat;
}"})

  # extract numeric allele counts for A:T:C:G
  syncCnts <- lapply(syncDt[,-1:-3,with=FALSE], function(col) {
    Sync2Cnts(col)
  } )

  # get rid of data that is no longer needed
  chr <- as.character(syncDt$V1)
  pos <- syncDt$V2
  ref <- syncDt$V3
  popCnt <- ncol(syncDt)-3
  rm(syncDt)
  gc()

  # sum up counts across populations (time points and replicates)
  sumCnts <- Reduce("+", syncCnts)

  # if only biallelic sites should be kept then identify those
  # else identify polymorphic sites in general
  ikeep = if(keepOnlyBiallelic) which(rowSums(sumCnts > 0) == 2) else which(rowSums(sumCnts > 0) > 1)
  # remove data from all other sites
  sumCnts <- sumCnts[ikeep, , drop=FALSE]
  chr <- chr[ikeep]
  pos <- pos[ikeep]
  ref <- ref[ikeep]

  # determine allele ranks for each position
  alleleRank <- rowRanks(sumCnts+rep(1:4/5,each=nrow(sumCnts)))
  rm(sumCnts)
  gc()
  # extract 2 most common alleles (considering all populations)
  alleleCnts <- lapply(syncCnts, function(pop) {
    cbind(major=t(pop[ikeep, ])[t(alleleRank == 4)],
          minor=t(pop[ikeep, ])[t(alleleRank == 3)])
  } )
  rm(syncCnts)
  gc()

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
  rm(alleleRank)
  gc()

  # if polarization is 'reference' -> check if ref-allele is either major or minor -> warning if not and polarization for minor instead
  if(polarization == "reference" && any(alleles$ref != alleles$minor & alleles$ref != alleles$major)) {
    warning("Cannot polarize for reference allele, because it is not among the two most common alleles for some SNPs. Changing polarization to 'minor'.")
    polarization <- "minor"
  }

  # combine generation and replicate info
  popInfo <- data.table(pop=1:popCnt, gen=gen, repl=repl)

  # add allele frequency and sequence coverage for each population
  for(r in unique(repl)) {
    for(i in seq(1:nrow(popInfo))[popInfo$repl == r]) {
      seqCov <- rowSums(alleleCnts[[i]])
      # compute allele frequencies according to 'polarization'
      if(polarization == "minor" || polarization == "rising") {
        alleles[,paste0("F", popInfo$gen[i], ".R", r, ".freq"):=alleleCnts[[i]][,"minor"]/seqCov]
      } else {
        alleles[,paste0("F", popInfo$gen[i], ".R", r, ".freq"):=ifelse(minor == ref, alleleCnts[[i]][,"minor"]/seqCov, alleleCnts[[i]][,"major"]/seqCov)]
      }
      alleles[,paste0("F", popInfo$gen[i], ".R", r, ".cov"):=seqCov]
    }
  }
  rm(alleleCnts)
  gc()

  # if required then polarize allele counts for the rising allele
  if(polarization == "rising" && length(unique(popInfo$gen)) > 1) {
    minGen <- min(popInfo$gen)
    maxGen <- max(popInfo$gen)

    # compute AFC between mean AF at first and last generation
    meanAF <- rowMeans(alleles[,grepl(paste0("F", maxGen, "\\.R[0-9]+\\.freq"), colnames(alleles)),with=FALSE], na.rm=TRUE) -
      rowMeans(alleles[,grepl(paste0("F", minGen, "\\.R[0-9]+\\.freq"), colnames(alleles)),with=FALSE], na.rm=TRUE)

    # polarize allele frequencies
    changePolarization <- isTRUE(meanAF < 0)
    for(pop in grep(".freq", colnames(alleles), value=TRUE, fixed=TRUE)) {
      alleles[[pop]][changePolarization] <- 1-alleles[[pop]][changePolarization]
    }

    # set column with rising allele
    alleles$rising <- ifelse(changePolarization, alleles$major, alleles$minor)
    rm(meanAF, changePolarization)
    gc()
  }

  # return sync-object for loaded sync-file
  return(new(Class="sync",
             gen=as.numeric(sub("F([0-9]+)\\.R[0-9]+.*", "\\1", colnames(alleles)[-1:-6])),
             repl=as.numeric(sub(".*\\.R([0-9]+)\\..*", "\\1", colnames(alleles)[-1:-6])),
             isAF=grepl(".*\\.freq$", colnames(alleles)[-1:-6]),
             polarization=polarization,
             alleles=alleles))
  }
