#'@importFrom MASS ginv
#'@importFrom stats median
#'@importFrom stats na.omit
#'@importFrom stats var
#'@import dplyr
#'@import readr
#'@import tidyverse

pairwise_dis<-function(assigndata,p,K,method){

  dis<-dis_K<-NULL
  for (s in 1:K) {
    for (t in 1:K) {
      if (t>s)
      {
        pairwise<-assigndata[which(!is.na(assigndata$assignment)),]

        if (p>1){
          u<-colSums((pairwise[which(pairwise$assignment==s),-1]))-
            colSums(as.matrix(pairwise[which(pairwise$assignment==t),-1]))
          cov<-cov(as.matrix(pairwise[,-1]))}

        else{
          u<-sum((pairwise[which(pairwise$assignment==s),-1]))-
            sum(as.matrix(pairwise[which(pairwise$assignment==t),-1]))
          cov<-var(as.matrix(pairwise[,-1]))}

        dis<-t(u)%*%ginv(cov)%*%u
        dis_K<-c(dis_K,dis)

      }
    }
  }

  if (method=='mean'){
    all_dis<-mean(dis_K)
  }
  if(method=='max'){
    all_dis<-max(dis_K)
  }
  if(method=='median'){
    all_dis<-median(dis_K)
  }
  if(method=='none'){
    all_dis<-dis_K
  }
  return(all_dis)
}

circle_random<-function(assigndata,start,d,K,p,q,method,n){
  for (j in start:(n-1)){
    num<-data.frame(matrix(NA,1,K))
    for (m in 1:K){
      num[1,m]<-nrow(assigndata[which(assigndata$assignment==m),])
    }
    if ((max(num)-min(num))<d)
    {ma<-data.frame(matrix(NA,1,K))
    for (k in 1:K) {
      assigndata$assignment[j+1]<-k
      ma[1,k]<-pairwise_dis(assigndata,p,K,method)
    }
    assign<-as.numeric(parse_number(colnames(ma[which.min(ma)])))
    all<-c(seq(1,K,1))
    assigndata$assignment[j+1]<-sample(c(assign,setdiff(all,assign)),
                                prob=c(q,rep((1-q)/(K-1),(K-1))),1,replace = TRUE)
    }
    else
    { len<-length(which(num==min(num)))
      min_vector<-which(num==min(num))
      if (len==1){
        assigndata$assignment[j+1]<-min_vector
      }
      else{
        assigndata$assignment[j+1]<-sample(x=c(min_vector),
                                           prob = c(rep(1/len,len)),
                                           1,replace = TRUE)
      }
    }
  }
  return(assigndata)
}


standard_pairwise_dis<-function(assigndata,p,K,method){
  dis<-dis_K<-NULL
  for (s in 1:K) {
    for (t in 1:K) {
      if (t>s)
      {

        pairwise<-assigndata[which(!is.na(assigndata$assignment)),]

        if (p>1){
          u<-colMeans((pairwise[which(pairwise$assignment==s),-1]))-
            colMeans(as.matrix(pairwise[which(pairwise$assignment==t),-1]))
          cov<-cov(as.matrix(pairwise[,-1]))}

        else{
          u<-mean((pairwise[which(pairwise$assignment==s),-1]))-
            mean(as.matrix(pairwise[which(pairwise$assignment==t),-1]))
          cov<-var(as.matrix(pairwise[,-1]))}

        dis<-2/K/K*nrow(pairwise)*(t(u)%*%ginv(cov)%*%u)
        dis_K<-c(dis_K,dis)


      }
    }
  }

  if (method=='mean'){
    all_dis<-mean(dis_K)
  }
  if(method=='max'){
    all_dis<-max(dis_K)
  }
  if(method=='median'){
    all_dis<-median(dis_K)
  }
  if(method=='none'){
    all_dis<-dis_K
  }
  return(all_dis)
}

circle_random_scr<-function(assigndata,start,d,K,p,q1,q2,method,n){
  for (j in start:(n-1)){
    num<-data.frame(matrix(NA,1,K))
    for (m in 1:K){
      num[1,m]<-nrow(assigndata[which(assigndata$assignment==m),])
    }
    if ((max(num)-min(num))<d)
    {ma<-data.frame(matrix(NA,1,K))
    for (k in 1:K) {
      assigndata$assignment[j+1]<-k
      ma[1,k]<-pairwise_dis(assigndata,p,K,method)
    }
    assign<-as.numeric(parse_number(colnames(ma[which.min(ma)])))
    all<-c(seq(1,K,1))
    assigndata$assignment[j+1]<-sample(c(assign,setdiff(all,assign)),
                                       prob=c(q1,rep((1-q1)/(K-1),(K-1))),1,replace = TRUE)
    }
    else
    { len<-length(which(num==min(num)))
    min_vector<-which(num==min(num))
    if (len==1){
      assigndata$assignment[j+1]<-sample(x=c(min_vector,setdiff(all,min_vector)),
                                        prob = c(q2,1-q2),
                                       1,replace = TRUE)
    }
    #else{
      #assigndata$assignment[j+1]<-sample(x=c(min_vector),
    #                                     prob = c(rep(1/len,len)),
     #                                    1,replace = TRUE)
    #}
    }
  }
  return(assigndata)
}


