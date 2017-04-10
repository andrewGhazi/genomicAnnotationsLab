#This code was originally written by Ian Campbell

#' Score a set of training genes based on an annotation source.
#' @param AnnotationMatrix a matrix with genes as rows and annotation variables as columns (with appropriate row/column names)
#' @param TrainingGenes a character vector of genes
#' @param MinPerAnno The minimum require number of genes annotated for an annotation to be considered
#' @param OR minimum Odds Ratio to be reported
#' @param PVal minimum p-value to be reported
#' @param MinAnno minimum number of annotations to select
#' @param MaxAnno maximum number of annotations to select
#'
#' @return a list with two elements: a named vector of Scores and a matrix of EnrichedAnnotations
#' @import Matrix
#' @importFrom corpcor fast.svd
#' @export
makeScore<-function(AnnotationMatrix=NULL,TrainingGenes=NULL,MinPerAnno=5, OR=1.5,PVal=.01,MinAnno=5,MaxAnno=200)
{
  aset <- TrainingGenes
  A <- AnnotationMatrix
  ### Compute Annotations Enriched in Gene Set
  enrich<-makeEnrichment(aset,
                         adj = A)

  ### Identify Annotations Fullfilling Criteria
  ###EDITS DONE HERE BY ANDREW GHAZI NOV 11 2016 # fixing the references to the output of makeEnrichment
  okcat<-which(as.numeric(as.matrix(enrich[,'Count']))>=MinPerAnno & as.numeric(as.matrix(enrich[,'OddsRatio']))>=OR & as.numeric(as.matrix(enrich[,'OddsRatio']))!=Inf & as.numeric(as.matrix(enrich[,'Pvalue']))<=PVal)

  ### Identify Most Significan Annotations
  if (length(okcat)>MaxAnno)
  {
    best<-okcat[order(-1*as.numeric(as.matrix(enrich[okcat,4])))][1:MaxAnno]
    okcat<-best
  }

  print(paste('Done selecting Annotations:',length(okcat)))

  if (length(okcat)<MinAnno) return('ERROR: Too few categories selected as enriched. Refine paremeters or increase the gene set.')

  ### Identify Training Set in Annotation Matrix
  aset.ind<-match(aset,rownames(A))
  if (sum(is.na(aset.ind))>0) aset.ind<-aset.ind[!is.na(aset.ind)]

  ### Determine Frequency Among Non-Training Genes
  bg.fq<-apply(A[-aset.ind,][,okcat],2,mean)
  print(paste('Done computing BG scores.'))

  ### Determine Frequency Among Training Genes
  aset.fq<-apply(A[aset.ind,][,okcat],2,mean)

  ### Calculate Background Inverse Covariance Matrix
  sub.cov=SimplestCmat2(A[,okcat])
  cov.mat<-fast.svd(sub.cov)
  acut<-1e-6
  passes.pos<-which(cov.mat$d>=acut)
  acovinv.bg<-cov.mat$u[,passes.pos]%*%diag(1/cov.mat$d[passes.pos])%*%t(cov.mat$u[,passes.pos])
  all.pos.score<-sapply(1:nrow(A), computeQform, aadj=A[,okcat], acovi=acovinv.bg, abg=bg.fq)
  print(paste('Done computing background quadratic form.'))

  ### Calculate Training Gene Inverse Covariance Matrix
  sub.cov=SimplestCmat2(A[aset.ind,][,okcat])
  cov.mat<-fast.svd(sub.cov)
  acut<-1e-6
  passes.pos<-which(cov.mat$d>=acut)
  acovinv.aset<-cov.mat$u[,passes.pos]%*%diag(1/cov.mat$d[passes.pos])%*%t(cov.mat$u[,passes.pos])
  set.pos.score<-sapply(seq(1:nrow(A)),computeQform,aadj=A[,okcat],acovi=acovinv.aset,abg=aset.fq)
  print(paste('Done computing training set quadratic form.'))

  ### Calculate Score Ratio
  R=all.pos.score/set.pos.score
  arange<-range(R,na.rm=T)
  Rset<-R[aset.ind]
  Rbg<-R[-aset.ind]
  if (min(R[aset.ind])> arange[1]) Rset<-c(R[aset.ind],arange[1])
  if (min(R[-aset.ind])>arange[1]) Rbg<-c(R[-aset.ind],arange[1])
  if (max(R[aset.ind])< arange[2]) Rset<-c(Rset,arange[2])
  if (max(R[-aset.ind])<arange[2]) Rbg<-c(Rbg,arange[2])

  ### Compute Tail Probability Ratios
  tail.set<-approx(sort(Rset),(1-(ecdf(Rset)(sort(Rset)))),xout=R)$y
  tail.bg<-approx(sort(Rbg),(1-(ecdf(Rbg)(sort(Rbg)))),xout=R)$y
  Z<-log2(tail.set/tail.bg)
  Z[is.infinite(Z)]<-min(Z[!is.infinite(Z)])/2
  Z<-(Z-mean(Z,na.rm=T))/sd(Z,na.rm=T)
  names(Z) <- rownames(A)

  ### Scores and Return Enriched Categories
  inner.output = list(Scores=Z,Enriched.Annotations=enrich[okcat,])
}

###### Computes Quadratic Form from Adjacency Matrix and Inverse Covariance Matrix

computeQform<-function(aindex,aadj,acovi,abg)
{

  arow<-aadj[aindex,]
  (t(arow-abg) %*%acovi %*% (arow-abg))/1
}

#' Detect enriched categories in an annotation source
#' @param gsyms a set of training genes
#' @param adj an annotation matrix
#' @param nom an optional vector of term names corresponding to the column names of the annotation matrix
#' @import Matrix
#' @importFrom corpcor fast.svd
makeEnrichment <- function (gsyms, adj = NULL, nom = NULL)
{
  express <- rownames(adj)
  okrows <- match(express, rownames(adj))
  okrows <- okrows[!is.na(okrows)]
  hits <- match(gsyms, rownames(adj))
  hits <- hits[!is.na(hits)]
  gsyms <- gsyms[which(hits %in% okrows)]
  hits <- hits[which(hits %in% okrows)]

  # ANDREW GHAZI EDIT NOV 16 2016
  # total.counts <- as.numeric(unlist(lapply(seq(1:ncol(adj)),
  #                                          function(i) {
  #                                            sum(adj[okrows, i], na.rm = T)
  #                                          })))
  total.counts = colSums(adj)

  asub <- adj[hits, ]
  rn <- rownames(asub)
  gn <- apply(asub, 2, function(x) {
    paste(rn[which(x == 1)], collapse = ",")
  })

  # ANDREW GHAZI EDIT NOV 16 2016
  #obs.counts <- as.numeric(as.matrix(apply(asub, 2, sum)))
  obs.counts = colSums(asub)

  expect <- length(hits) * total.counts/length(okrows)
  N <- length(okrows)
  an <- length(gsyms)
  atest <- t(lapply(1:length(obs.counts), function(x) {
    ans <- try(fisher.test(rbind(c(obs.counts[x], an - obs.counts[x]),
                                 c(total.counts[x] - obs.counts[x], N - total.counts[x] -
                                     an + obs.counts[x]))), silent = TRUE)
    if (!(class(ans) == "try-error")) {
      return(c(ans$estimate, ans$p.value))
    }
    else {
      return(c(NA, NA))
    }
  }))
  atest <- t(sapply(atest, function(x) {
    x[1:2]
  }))
  fin <- cbind(colnames(adj), obs.counts, total.counts, (atest),
               gn)

  if (!is.null(nom)) {
    try(nmap <- nom[match(fin[, 1], nom[, 1]), 2], silent = T)
    if (class(nmap) == "try-error")
      stop("ERROR IN NAMES")
    fin <- cbind(fin[, 1], nmap, fin[, -1])
    colnames(fin) <- c("TermID", "TermName", "Count", "Total",
                       "OddsRatio", "Pvalue", "Genes")
  } else {
    colnames(fin) <- c("TermID", "Count", "Total",
                       "OddsRatio", "Pvalue", "Genes")
  }


  return(fin)
}

SimplestCmat2<-function(aadj,n=1)
{
  m<-nrow(aadj)
  alloverlap<-t(aadj)%*%aadj
  maincount<-diag(alloverlap)
  pairproduct<-maincount%*%t(maincount)
  leading.coefficient<--n/m^2*(m-n)/(m-1)
  cov.mat<-leading.coefficient * (pairproduct - m*alloverlap)
  vars<- maincount*(m-maincount) * (-1*leading.coefficient)
  diag(cov.mat)<-vars
  rownames(cov.mat)<-colnames(cov.mat)<-colnames(aadj)
  cov.mat
}

makeCrossValidScore <- function(AnnotationMatrix = NULL,TrainingGenes = NULL,LeaveOneOut = FALSE, Fold = 10, Multi.Core = FALSE, N.Cores = 4, ...){
  MyAnnoMat <- AnnotationMatrix
  ###Split Training Genes
  if(LeaveOneOut == TRUE){
    Gene.List <- split(TrainingGenes,seq(1,length(TrainingGenes)))
  }
  else {
    Gene.List <- split(TrainingGenes,floor(seq_along(TrainingGenes)/ceiling(length(TrainingGenes)/Fold)))
  }

  if(Multi.Core == TRUE){
    require(parallel)
    Scores <- mclapply(Gene.List,function(x){
      ### Remove Cross Validation Set
      Gene.List <- unlist(Gene.List)
      NewTraining <- Gene.List[!Gene.List %in% x]

      ### Calculate Scores
      Result <- makeScore(AnnotationMatrix=MyAnnoMat,TrainingGenes=NewTraining,...)$Scores

      ###Retreive Gene Scores
      Gene.Scores <- Result[match(x,names(Result))]
      Gene.Percentiles <- ecdf(Result)(Gene.Scores)
      names(Gene.Percentiles) <- x
      return(Gene.Percentiles)
    },mc.cores=N.Cores)
  }
  else{
    Scores <- lapply(Gene.List,function(x){
      ### Remove Cross Validation Set
      Gene.List <- unlist(Gene.List)
      NewTraining <- Gene.List[!Gene.List %in% x]

      ### Calculate Scores
      Result <- makeScore(AnnotationMatrix=MyAnnoMat,TrainingGenes=NewTraining,...)$Scores

      ###Retreive Gene Scores
      Gene.Scores <- Result[match(x,names(Result))]
      Gene.Percentiles <- ecdf(Result)(Gene.Scores)
      names(Gene.Percentiles) <- x
      return(Gene.Percentiles)
    })
  }

  Scores <- unlist(Scores)
  try(names(Scores) <- TrainingGenes)
  CV.Result <- list()
  class(CV.Result) <- "GeneAnnoCV"
  if(LeaveOneOut == TRUE){
    CV.Result$Type <- "Leave One Out"
  }
  else{
    CV.Result$Type <- paste(Fold,"Fold")
  }
  CV.Result$Scores <- Scores
  return(CV.Result)
}

print.GeneAnnoCV <- function(x, ...){
  cat("$Type\n  Type of Cross Validation:\n" , x$Type , "\n\n")
  cat("$Scores\n  Cross Validated Scores:\n" , x$Scores , "\n")
}

plot.GeneAnnoCV <- function(x,AUC=TRUE, ...){
  TrainingPer <- x$Scores
  Xs <- seq(0,1,by=.005)
  Ys <- sapply(Xs,function(x){sum(TrainingPer > 1-x,na.rm=T)/length(TrainingPer[!is.na(TrainingPer)])})
  plot(Xs,Ys,type="l",xlab="1 - Score Percentile", ylab="Fraction of Training Genes", ...)
  abline(0,1)
  Approx.AUC <- sum(diff(Xs)*(Ys[-1]+Ys[-length(Ys)]))/2
  try(if(AUC){text(0.75,0.1,paste("AUC:",round(Approx.AUC,3)))},silent=T)
}

lines.GeneAnnoCV <- function(x, ...){
  TrainingPer <- x$Scores
  Xs <- seq(0,1,by=.005)
  Ys <- sapply(Xs,function(x){sum(TrainingPer > 1-x,na.rm=T)/length(TrainingPer[!is.na(TrainingPer)])})
  lines(Xs,Ys, ...)
}
