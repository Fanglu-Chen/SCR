#'@export
#'@title Implement of SMART
#'
#'@description randomize patients into treatment groups by SMART
#'@param covariate A dataframe of patients' covariates which are needed to balance, the first column is patients' unique ID number and the rest of the columns are corresponding patients' covariates.
#'@param assignment A dataframe of patients' allocation which has two columns; If you had allocated partial patients, please input them (the first column is patients' unique ID number and the second colunmn is patients' allocation), please pre-set the allocation group number and enter the number in the second column. IF all the patients are not be allocated, please input 'assignment = NA' directly.
#'@param K An integer; number of assignments of the trial.
#'@param d An integer; restriction for sample sizes.
#'@param q A number; allocation probability, default = 0.75
#'@param method If K=2, please IGNORE this parameter and don't input anything; If K>2, you can input one of these texts: 'mean', 'max' or 'median'.

#'@examples
#'#simulate covariates of patients
#'p=6;n=30
#'sigma<-diag(p);mean<-c(rep(0,p))
#'data <- mvrnorm(n, mean, sigma)
#'covariate<-as.data.frame(data)
#'#IF all the patients are not be allocated
#'smart(covariate = covariate,assignment = NA,K=3,d=5,q=0.75)
#'#IF you had allocated partial patients
#'assignment<-c(1,2,2)
#'smart(covariate = covariate,assignment = assignment,K=3,d=3,q=0.75)

smart<-function(covariate, assignment, K, d, q=0.75,method='mean'){
  n<-nrow(covariate)
  p<-ncol(covariate)

  if (is.na(assignment)[1]){
    assignment<-data.frame(assignment=rep(NA,n))
  }

  else{
    aln<-length(assignment)
    assignment<-data.frame(assignment=c(assignment,rep(NA,n-aln)))
  }

  assigndata<-cbind(data.frame(assignment),
                    data.frame(covariate))
  names(assigndata)[1]<-'assignment'
  assigndata<-dplyr::arrange(assigndata, is.na(assigndata$assignment))


  if (nrow(assigndata[is.na(assigndata$assignment),])==n)
  {
    assigndata$assignment<-c(seq(1,K,1),rep(NA,(n-K)))
    start<-K
    assigndata<-circle_random(assigndata,start,d,K,p,q,method,n)
  }
  else
  {
    noassign<-setdiff(seq(1,K,1),unique(na.omit(assigndata$assignment)))
    if (length(noassign)==0){
      start<-nrow(assigndata[!is.na(assigndata$assignment),])
      assigndata<-circle_random(assigndata,start,d,K,p,q,method,n)
    }
    else
    {
      assigndata$assignment<-c(assigndata[!is.na(assigndata$assignment),'assignment'],noassign,
                               rep(NA,(n-length(assigndata[!is.na(assigndata$assignment),'assignment'])-
                                         length(noassign))))
      start<-nrow(assigndata[!is.na(assigndata$assignment),])
      assigndata<-circle_random(assigndata,start,d,K,p,q,method,n)
    }
  }
  R = NULL
  R$assignment<-assigndata[,1]
  R$sample_size<-as.data.frame(assigndata%>%group_by(assignment)%>%count(assignment))
  R$Mahalanobis_Distance<-standard_pairwise_dis(assigndata,p,K,method)
  return(R)
  }
