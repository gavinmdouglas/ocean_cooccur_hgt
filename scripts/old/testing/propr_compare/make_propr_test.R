suppressMessages(library(data.table))
suppressMessages(library(compositions))
suppressMessages(library(dplyr))
suppressMessages(library(propr))
suppressMessages(library(zCompositions))
suppressMessages(library(parallel))
suppressMessages(library(doParallel))

my.ivar2index <- function (counts, ivar)
{
  if (missing(ivar))
    ivar <- 0
  if (!is.vector(ivar))
    stop("Provide 'ivar' as vector.")
  `%is%` <- function(a, b) identical(a, b)
  if (ivar %is% 0 | ivar %is% NA | ivar %is% NULL | ivar %is%
      "all" | ivar %is% "clr" | ivar %is% "tss") {
    use <- 1:ncol(counts)
  }
  else if (ivar %is% "iqlr") {
    if (any(counts == 0)) {
      message("Alert: Replacing 0s with next smallest value.")
      zeros <- counts == 0
      counts[zeros] <- min(counts[!zeros])
    }
    counts.clr <- apply(log(counts), 1, function(x) {
      x - mean(x)
    })
    counts.var <- apply(counts.clr, 1, var)
    quart <- stats::quantile(counts.var)
    use <- which(counts.var < quart[4] & counts.var > quart[2])
  }
  else {
    if (is.character(ivar)) {
      if (!all(ivar %in% colnames(counts)))
        stop("Some 'ivar' not in data.")
      use <- which(colnames(counts) %in% ivar)
    }
    else {
      use <- sort(ivar)
    }
  }
  return(use)
}



my.propr <- function (counts,
                      min_occ = NA,
                      metric = "rho",
                      ivar = c("clr","tss","none"), select, symmetrize = FALSE,
                      pseudocount="none", p = 100, only_positive = FALSE) {

  # TSS relative abundances
  if (pseudocount == "none" & ivar == "tss" ) {
    metric = "cor"
  }
  message("* Pseudocount : ",pseudocount)
  message("* Log transfo : ",ivar)
  message("* Metric : ",metric)

  if (any(counts < 0) & ivar!="none") {
    stop("Data may not contain negative measurements.")
  }
  if (any(is.na(counts)))
    stop("Remove NAs from 'counts' before proceeding.")
  if ("data.frame" %in% class(counts))
    counts <- as.matrix(counts)
  if (is.null(colnames(counts)))
    colnames(counts) <- as.character(1:ncol(counts))
  if (is.null(rownames(counts)))
    rownames(counts) <- as.character(1:nrow(counts))
  ct <- counts

  # pseudocount
  if(ivar=="clr" && any(as.matrix(counts) == 0)) {
    # message(paste("Alert: Replacing 0s with ",pseudocount,".",sep=""))
    if (pseudocount == "bayesian") {
      mat.ct <- data.matrix(ct)
      ct <- as.matrix(zCompositions::cmultRepl(X = mat.ct, output="prop"))
      colnames(ct) <- as.character(1:ncol(counts))
      rownames(ct) <- as.character(1:nrow(counts))
    } else if (pseudocount == "multiplicative_replacement") {
      zeros <- ct == 0
      pseudocount = min(ct[ct!=0]) * 0.65
      ct[zeros] <- pseudocount
    } else if (pseudocount == "add_one_everywhere") {
      ct = ct + 1
    } else {
      stop("Unknown pseudocount !")
    }
  }
  print(paste(dim(ct)[2],"species and",dim(ct)[1],"samples"))

  # transformation
  if (ivar == "clr") {
    use <- my.ivar2index(ct, ivar)
    # message("Alert: Saving log-ratio transformed counts to @logratio.")
    logX <- log(ct)
    logSet <- logX[, use, drop = FALSE]
    ref <- rowMeans(logSet)
    lr <- sweep(logX, 1, ref, "-")
    if (min_occ != 0) {
      ct_bool=counts
      ct_bool[ct_bool!=0]=1
      sp_filter = which(apply(ct_bool,2,sum)>=min_occ)
      message(length(sp_filter)," species remaining out of ",
              dim(ct)[1]," with min occ = ",min_occ)
      # print(sp_filter)
      print(dim(ct))
      lr = lr[,sp_filter]
      ct = ct[,sp_filter]
      print(dim(ct))
      #message("dimensions of ct: ",dim(ct)[1]," samples and ",dim(ct)[2]," species")
    }
    if (only_positive) { lr <- lr + abs(min(lr)) ; message("Alert: Only positive CLR values")}

  } else if (ivar=="tss") {
    # message("Alert: Skipping the log-ratio transformation.")
    # alpha <- NA
    ct <- t(apply(counts,1,function(i) i/sum(i)))
    ct[is.na(ct)] = 0
    lr <- ct
  } else { # no transfo
    lr <- ct
  }

  lrv <- lr2vlr(lr)


  # rho or cor
  metric <- metric[1]
  if (metric == "rho") {
    ## ajout d'une valeur fixe
    # mini = abs(min(lr)); lr = lr + mini
    mat <- lr2rho(lr)
  } else if (metric == "cor"){
    mat <- stats::cor(lr)
  } else {
    stop("Provided 'metric' not recognized.")
  }

  # Build propr object
  result <- new("propr")

  result@alpha <- as.numeric(NA)
  result@metric <- metric #rho
  result@ivar <- ivar #clr
  # print(paste("min lr :",min(lr)))
  result@logratio <- as.data.frame(lr)
  result@pairs <- vector("numeric")

  result@matrix <- mat
  colnames(result@matrix) <- colnames(result@logratio)
  rownames(result@matrix) <- colnames(result@logratio)

  # Set up @permutes
  result@permutes <- list(NULL)
  if(p > 0){
    # Shuffle row-wise so that per-sample CLR does not change
    # message("Alert: Fixing permutations to active random seed.")
    permutes <- vector("list", p)
    for(ins in 1:length(permutes)) permutes[[ins]] <- t(apply(ct, 1, sample))
    result@permutes <- permutes
  }

  # Set up @results
  labels <- labRcpp(ncol(lr))
  result@results <-
    data.frame(
      "Partner" = labels[[1]],
      "Pair" = labels[[2]],
      "lrv" = lltRcpp(lrv),
      "metric" = factor(metric),
      "alpha" = factor(NA),
      "propr" = lltRcpp(mat)
    )

  # Initialize @results -- Tally frequency of 0 counts
  if (min_occ!=0) {
    counts=counts[,sp_filter]
  }
  if(any(as.matrix(counts) == 0) & !is.na(ivar)){
    # message("Alert: Tabulating the presence of 0 counts.")
    result@results$Zeros <- ctzRcpp(as.matrix(counts)) # count 0s
  }

  # message("Alert: Use '[' to index proportionality matrix.")
  # message("Alert: Use 'updateCutoffs' to calculate FDR.")

  return(result)
}

my.updateCutoffs.propr <- function (object, min_occ=NA, cutoff, pseudocount, ncores)
{
  getFdrRandcounts <- function(ct.k) {
    pr.k <- suppressMessages(my.propr(ct.k, min_occ=min_occ,
                                      object@metric, pseudocount = "multiplicative_replacement",
                                      ivar = object@ivar))#, alpha = object@alpha))
    pkt <- pr.k@results$propr
    sapply(FDR$cutoff, function(cut) {
      if (object@metric == "rho" | object@metric == "cor") {
        sum(pkt > cut)
      }
      else {
        sum(pkt < cut)
      }
    })
  }
  if (identical(object@permutes, list(NULL)))
    stop("Permutation testing is disabled.")
  if (identical(cutoff, NA))
    return(object)
  FDR <- as.data.frame(matrix(0, nrow = length(cutoff), ncol = 4))
  colnames(FDR) <- c("cutoff", "randcounts", "truecounts",
                     "FDR")
  FDR$cutoff <- cutoff

  # index of neg and pos pairs
  neg_i = which(object@results$propr < 0)
  pos_i = which(object@results$propr > 0)

  p <- length(object@permutes)
  # si on parallèlise
  if (ncores > 1) {
    # packageCheck("parallel")
    cl <- parallel::makeCluster(ncores)
    parallel::clusterExport(cl=cl,list("my.propr","my.ivar2index"))
    #lapply(object@permutes,getFdrRandcounts)
    randcounts <- parallel::parLapply(cl = cl, X = object@permutes,
                                      fun = getFdrRandcounts)
    FDR$randcounts <- apply(as.data.frame(randcounts), 1,
                            sum)
    parallel::stopCluster(cl)
  }
  # si on ne parallèlise pas
  else {
    for (k in 1:p) {
      numTicks <- progress(k, p, numTicks)
      ct.k <- object@permutes[[k]]
      pr.k <- suppressMessages(my.propr(ct.k, metric = object@metric, min_occ=min_occ,
                                        ivar = object@ivar, pseudocount=pseudocount))#, alpha = object@alpha))
      pkt <- pr.k@results$propr
      for (cut in 1:nrow(FDR)) {
        if (object@metric == "rho" | object@metric ==
            "cor") {

          # neg or pos pairs ?
          if (FDR[cut, "cutoff"] < 0) {
            FDR[cut, "randcounts"] <- FDR[cut, "randcounts"] +
              sum(pkt[neg_i] < FDR[cut, "cutoff"])
          } else {
            FDR[cut, "randcounts"] <- FDR[cut, "randcounts"] +
              sum(pkt[pos_i] > FDR[cut, "cutoff"])
          }


        } # pour chaque cut off, on compte le nombre de propr qui sont strict. supérieurs dans l'objet permuté donc random !
        else {
          FDR[cut, "randcounts"] <- FDR[cut, "randcounts"] +
            sum(pkt < FDR[cut, "cutoff"]) # osef
        }
      }
    }
  }
  FDR$randcounts <- FDR$randcounts/p # on divise par le nombre de permutations, donc moyenne sur l'ensemble des objets permutés
  for (cut in 1:nrow(FDR)) {
    if (object@metric == "rho" | object@metric == "cor") {

      # neg or pos pairs ?
      if (FDR[cut, "cutoff"] < 0) {
        FDR[cut, "truecounts"] <- sum(object@results$propr[neg_i] < FDR[cut, "cutoff"])
      } else {
        FDR[cut, "truecounts"] <- sum(object@results$propr[pos_i] > FDR[cut, "cutoff"])
      }


    }
    else {
      FDR[cut, "truecounts"] <- sum(object@results$propr <
                                      FDR[cut, "cutoff"])
    }
    FDR[cut, "FDR"] <- FDR[cut, "randcounts"]/FDR[cut, "truecounts"] # et le FDR c'est randcounts/trucounts
    # save step by step
    write.csv(as.data.frame(FDR),
              paste("motus-profiles_pos-cutoff_fdr.csv",sep=""),
              quote=F)
  }
  object@fdr <- FDR
  return(object)
}

environment(my.updateCutoffs.propr) <- environment(propr)
environment(my.propr) <- environment(propr)
environment(my.ivar2index) <- environment(propr)


df=read.csv("/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/metaG_rpkm_freeliv.tsv.gz", sep="\t",row.names=1, header = TRUE)
rownames(df) <- NULL

ct.propr_raw = my.propr(df,
                    min_occ=10,
                    metric="rho",
                    ivar= "clr",
                    pseudocount= "multiplicative_replacement")

ct.propr.results = ct.propr_raw@results
ct.propr.results=data.frame(ct.propr.results)

ct.propr = my.updateCutoffs.propr(ct.propr_raw,
                                  min_occ=10,
                                  cutoff=seq(0.05,0.95,0.05),
                                  pseudocount="multiplicative replacement",
                                  ncores=64)

ct.propr.results=data.frame(ct.propr.results)



ct.propr.results$Partner_name <- colnames(df)[ct.propr.results$Partner]
ct.propr.results$Pair_name <- colnames(df)[ct.propr.results$Pair]

write.table(x = ct.propr.results, file = "~/tmp/propr_freeliv_earlier_ms_code.tsv",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")