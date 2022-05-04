#' @title BNP Testing Procedure for Paired Samples function
#'
#'
#' @description
#' Given two vectors of numerical values, this function returns the
#' result for Bayesian nonparametric hypothesis testing for paired
#' samples proposed by Pereira et al. (2020), performing an analytical and graphical comparison of the
#' marginal distributions of the data set.
#'
#'
#' @param x a numeric vector of data values taken prior to measurement.
#' @param y a numeric vector of data values taken post measurement.
#' @param n.mcm the number of simulations for the MCMC in the Gibbs sampling (suggested: 10000).
#'
#'
#' @return A list with three items. The first element (\code{sampling.parameters})
#'  is a list of the parameters estimated by Gibbs sampling, which are returned
#'  as data frames within each of the iterations. The second element
#'  (\code{posterior.probability.H1}) refers to the posterior probability
#'  for the alternative hypothesis, i.e. that differences between the marginal
#'  distributions occur. The third element (\code{standardized.data}) refers
#'  to the original standardized data set, which will be useful when a when
#'  applying other functions of the package.
#'
#'
#' @importFrom pracma pinv
#' @importFrom MASS mvrnorm
#' @importFrom mvtnorm dmvnorm
#' @importFrom truncnorm rtruncnorm
#' @importFrom truncnorm dtruncnorm
#' @import plotly
#' @import ggplot2
#' @import dplyr
#'
#'
#' @note For a proper execution of the function it is required
#' that the data vectors have the same length, otherwise the function
#' will return an error message.
#'
#'
#' @examples
#'
#' \dontrun{
#' x <- rnorm(30,3,2)
#' y <- rnorm(30,4,3)
#' BNP.test(x, y, n.mcm=10000)
#'
#' y <- matrix(runif(300), ncol = 300)
#' y <- apply(y, 2, function(i) if (i < 0.5) {
#'   y <- rmvnorm(1, mean = c(0,-3), sigma = matrix(c(1,0.8,0.8,1),nrow = 2,byrow = T))
#' }else{
#'   y <- rmvnorm(1, mean = c(0,3), sigma = matrix(c(1,0.8,0.8,1),nrow = 2,byrow = T))
#' })
#' y <- t(y)
#' x1 <- y[,1]
#' y1 <- y[,2]
#' BNP.test(x1, y1, n.mcm=10000)
#'}
#'
#' @references{
#' Pereira, L. A., Taylor-Rodriguez, D. & Gutierrez, L. (2020), \emph{A Bayesian nonparametric
#' testing procedure for paired samples}. Biometrics 76(1), 1-14.
#' }
#'
#'
#' @export
BNP.test <- function(x, y, n.mcm){

  # Data loading ---------------------------------------------------------------

  x.data<-x
  y.data<-y

  if(length(x.data)!=length(y.data)){
    stop("Error: the length of the vectors x, y is not equal")
  }

  sample.y<-cbind(x.data,y.data)


  # Definition of variables and initial values ------------------------------

  mu.0<-c(0,0)

  a2<-5
  b2<-1

  v0<-0.1

  epsilon<-0.01

  a3<-50
  b3<-5

  a0<-0.01
  b0<-0.01

  a.gamma<-(1/2)
  b.gamma<-(1/2)

  a1<-20
  b1<-1

  kappa<- 10

  param1.sigma2<-0.1*(1/v0)
  param2.sigma2<-0.1*(1/v0)

  param2.sigma2.alt<-0.1
  param1.sigma2.alt<-0.1

  s<- param1.sigma2.alt

  k<-1:2

  if(nrow(sample.y)%%2 == 0){
    cluster.0 <- rep(c(1,2),nrow(sample.y)/2)
  }else{
    cluster.0 <- rep(c(1,2),(nrow(sample.y)+1)/2)
    cluster.0 <- cluster.0[-(nrow(sample.y)+1)]
  }

  M<-5
  hypothesis<-0
  pi.0<-0.5

  z1<-rep(1,nrow(sample.y))
  z2<-rep(1,nrow(sample.y))
  delta<-rep(0,nrow(sample.y))

  data.initial<-data.frame(sample.y,cluster.0,z1,z2,delta)

  b<-matrix(sort(unique(cluster.0)),ncol = 1)
  parameters<-matrix(rep(0,length(b)*7),ncol=length(b))
  rownames(parameters)<-c("Beta1","Beta2","variance1","variance2","gamma1","gamma2","tau")

  parameters[1,]<-rep(mean(data.initial[,1]),length(b))
  parameters[2,]<-rep(mean(data.initial[,2])-mean(data.initial[,1]),length(b))
  parameters[3,]<-rep(var(data.initial[,1]),length(b))
  parameters[4,]<-rep(var(data.initial[,2]),length(b))
  parameters[5,]<-rep(0.01,length(b))
  parameters[6,]<-rep(0.02,length(b))
  parameters[7,]<-rep(3,length(b))


  # Gibbs sampling function for posterior inference of the atoms -------------

  atoms<-function(b,mat.data,t){
    z1.aux<-rep(0,nrow(mat.data))
    z2.aux<-rep(0,nrow(mat.data))
    vector.delta.i.aux<-rep(0,nrow(mat.data))

    y1<-mat.data[mat.data[,3]==b,1]
    y2<-mat.data[mat.data[,3]==b,2]
    z1<-mat.data[mat.data[,3]==b,4]
    z2<-mat.data[mat.data[,3]==b,5]
    sum.y1<-sum(y1)
    sum.y2<-sum(y2)
    n<-length(y1)

    Beta1<-t[1,b]
    Beta2<-t[2,b]
    variance1<-t[3,b]
    variance2<-t[4,b]
    gamma1<-t[5,b]
    gamma2<-t[6,b]
    tau<-t[7,b]

    X<-matrix(c(1,1,0,1),ncol=2)

    ## Updating the mean of conditional distributions

    sigma0<-diag(c(10,kappa),nrow=2)

    product.0<-solve(sigma0) %*% mu.0

    matrix_var<-pracma::pinv(n*t(X)%*%pracma::pinv(matrix(c(variance1 +tau ,rep((1-2*gamma1)*(1-2*gamma2)*tau,2), tau +(variance1*variance2)),byrow=F,ncol=2))%*%X + solve(sigma0))
    matrix_mu<-matrix_var %*% (product.0 + t(X)%*%pracma::pinv(matrix(c(variance1 +tau ,rep((1-2*gamma1)*(1-2*gamma2)*tau,2), tau +(variance1*variance2)),byrow=F,ncol=2))%*%matrix(c(sum.y1,sum.y2)))

    betas<-MASS::mvrnorm(n = 1,matrix_mu,matrix_var)
    Beta1<-betas[1]
    Beta2<-betas[2]

    ## Updating the covariance matrix of conditional distributions ##

    # Updating delta

    numerator<-(z1*(y1-Beta1)/variance1)+(z2*(y2-Beta1-Beta2)/(variance1*variance2))
    denominator<-(1/variance1)+(1/(variance2*variance1))+(1/tau)

    mean.delta<-numerator/denominator

    variance.delta<-1/denominator

    delta<-rnorm(n,mean.delta, sqrt(variance.delta))


    # Updating tau

    tau<-1/rgamma(1,shape=a0+(n/2), rate=(1/2)*(sum(delta^2)) + b0)


    # Updating sigma1

    variance1<-1/rgamma(1,shape=n+epsilon, rate=(1/2)*sum((y1-(Beta1+z1*delta))^2)+(1/(2*variance2))*sum((y2-(Beta1+Beta2+z2*delta))^2)+epsilon)


    # Updating sigma2

    if(hypothesis==1){

      param1.sigma2.alt<-s
      param2.sigma2.alt<-s

      param1.sigma2<-s/v0
      param2.sigma2<-s/v0

    }else{


      param1.sigma2.alt<-s*v0
      param2.sigma2.alt<-s*v0

      param1.sigma2<-s
      param2.sigma2<-s

    }
    if(hypothesis==1){
      param1.final.sigma2<-param1.sigma2.alt
      param2.final.sigma2<-param2.sigma2.alt
    }else{
      param1.final.sigma2<-param1.sigma2
      param2.final.sigma2<-param2.sigma2
    }

    alpha.sigma2<-param1.final.sigma2+(n/2)
    beta.sigma2<-(1/(2*variance1))*sum((y2-(Beta1+Beta2+z2*delta))^2) + param2.final.sigma2

    variance2<-1/rgamma(1,shape=alpha.sigma2, rate=beta.sigma2)


    # Updating z1

    prob.post.z1_neg<-gamma1*dnorm(y1,mean=Beta1-delta,sd=sqrt(variance1))
    prob.post.z1_positivo<- (1-gamma1)*dnorm(y1,mean=Beta1+delta,sd=sqrt(variance1))
    vec.prob_z1<-matrix(c(prob.post.z1_neg/(prob.post.z1_neg+prob.post.z1_positivo), prob.post.z1_positivo/(prob.post.z1_neg+prob.post.z1_positivo)),ncol=2)

    k<-matrix(c(1:n),ncol=1)
    z1<-apply(k,1, function(k) sample(c(-1,1),1,replace=F, prob=vec.prob_z1[k,]))


    # Updating z2

    prob.post.z2_neg<-gamma2*dnorm(y2,mean=Beta1+Beta2-delta,sd=sqrt(variance1*variance2))
    prob.post.z2_positivo<-(1-gamma2)*dnorm(y2,mean=Beta1+Beta2+delta,sd=sqrt(variance1*variance2))


    vec.prob_z2<-matrix(c(prob.post.z2_neg/(prob.post.z2_neg+prob.post.z2_positivo), prob.post.z2_positivo/(prob.post.z2_neg+prob.post.z2_positivo)),ncol=2)

    z2<-apply(k,1, function(k) sample(c(-1,1),1,replace=F, prob=vec.prob_z2[k,]))


    # Updating gamma.j

    gamma1<-rbeta(1, shape1=a.gamma+ sum(z1==-1), shape2=b.gamma-sum(z1==-1)+n)

    gamma2<-rbeta(1, shape1=a.gamma+ sum(z2==-1), shape2=b.gamma-sum(z2==-1)+n)


    z1.aux[mat.data[,3]==b]<-z1
    z2.aux[mat.data[,3]==b]<-z2
    vector.delta.i.aux[mat.data[,3]==b]<-t(delta)

    results<-c(Beta1,Beta2,variance1,variance2,gamma1,gamma2,tau)
    list(results,z1.aux,z2.aux,vector.delta.i.aux)
  }


  vec.zetas<-rep(NA,n.mcm)
  final.parameters<-list()
  densities_1<-list()
  densities_2<-list()
  min_value<-min(sample.y)-1
  max_value<-max(sample.y)+1
  y11<-seq(min_value,max_value,length.out = 200)
  y22<-seq(min_value,max_value,length.out = 200)


  # Simulation with MCMC algorithm ------------------------------------------

  for(a in 1:n.mcm){

    b<-matrix(sort(unique(cluster.0)),ncol = 1)
    theta1.j<-apply(b,1,atoms,mat.data=data.initial,t=parameters)

    for(i in 1:length(theta1.j)){
      delta<-delta+theta1.j[[i]][[4]]
    }

    z1.aux<-rep(0,nrow(data.initial))
    for(i in 1:length(theta1.j)){
      z1.aux<-z1.aux+theta1.j[[i]][[2]]
    }

    z2.aux<-rep(0,nrow(data.initial))
    for(i in 1:length(theta1.j)){
      z2.aux<-z2.aux+theta1.j[[i]][[3]]
    }

    data.initial[,4]<-z1.aux
    data.initial[,5]<-z2.aux
    data.initial[,6]<-delta

    parameters<-matrix(rep(0,length(b)*7),ncol=length(b))
    rownames(parameters)<-c("Beta1","Beta2","variance1","variance2","gamma1","gamma2","tau")

    for(i in 1:length(theta1.j)){
      parameters[,i]<-theta1.j[[i]][[1]]
    }


    # Updating the weights by stick-breaking process

    v.j<-apply(b,1,function(b) rbeta(1,1+sum(cluster.0==b),M + sum(cluster.0>b)))
    prod.cum<-cumprod(1-v.j)
    w<-v.j*c(1,prod.cum[-length(prod.cum)])


    if(sum(w)<1){
      p.j<-c(w,1-sum(w))
    }else{
      p.j<-w
    }

    if(ncol(parameters)<length(p.j)){


      sigma0<-diag(c(10,kappa),nrow=2)

      param1.final.sigma2<-s
      param2.final.sigma2<-s

      mu1<-MASS::mvrnorm(n = 1,mu.0,sigma0)
      variance1.j<-1/rgamma(1,shape=epsilon, rate=epsilon)
      variance2.j<-1/rgamma(1,shape=param1.final.sigma2, rate=param2.final.sigma2)
      tau.j<-1/rgamma(1,shape=a0, rate=b0)



      gamma1.j<-rbeta(1, shape1=a.gamma, shape2=b.gamma)
      gamma2.j<-rbeta(1, shape1=a.gamma, shape2=b.gamma)


      parameters<-cbind(parameters,c(mu1,variance1.j,variance2.j,gamma1.j,
                                     gamma2.j,tau.j))
    }


    # Updating the Mass parameter

    eta<-rbeta(1,M+1,nrow(data.initial))

    tau.eta<-(a1 + ncol(parameters)-1)/(nrow(data.initial)*b1 - nrow(data.initial)*log(eta)+a1+ncol(parameters)-1)


    u.M<-runif(1)
    if(u.M<tau.eta){
      M<-rgamma(1,shape=a1+ncol(parameters), rate=b1 - log(eta))
    }else{
      M<-rgamma(1,shape=a1+ncol(parameters)-1, rate=b1 - log(eta))
    }

    # Sampling the latent variable

    u<-rep(NA,nrow(data.initial))
    data1<-cbind(nrow(data.initial),cluster.0,u)
    g<-sort(unique(cluster.0))
    sort.clust<-matrix(g,ncol=1)

    for (i in 1:length(sort.clust)){
      data1[data1[,2]==sort.clust[i],3]<-runif(sum(cluster.0==sort.clust[i]),0,p.j[i])
    }

    u.j<-data1[,3]


    # Updating the bivariate density

    x<-matrix(c(1:length(p.j)),ncol=1)
    l<-apply(x,1,function (i) p.j[i] > u.j)*1
    l1<-apply(x,1,function (i) l[,i]*i)
    N<-apply(l1,1,max)
    k<-1:max(N)

    r<-matrix(k,ncol=1)
    indicadoras<-apply(r,1,function(j) (p.j[j]>u.j)*1)

    X<-matrix(c(1,1,0,1),ncol=2)

    densities<-apply(r,1, function(j) mvtnorm::dmvnorm(sample.y, mean=t(X%*%parameters[1:2,j]),
                                              sigma=matrix(c(parameters[3,j]+parameters[7,j],rep((1-2*parameters[5,j])*(1-2*parameters[6,j])*parameters[7,j],2),
                                                             (parameters[3,j]*parameters[4,j])+parameters[7,j]),ncol=2)))

    densities.1<-apply(matrix(1:length(p.j), ncol = 1),1,function(j) p.j[j] * dnorm(y11,parameters[1,j], sqrt(parameters[3,j]+parameters[7,j])))
    densities.2<-apply(matrix(1:length(p.j), ncol = 1),1,function(j) p.j[j]*dnorm(y22, parameters[1,j]+parameters[2,j],sqrt((parameters[3,j]*parameters[4,j])+parameters[7,j])))
    density.mixture1<-rowSums(densities.1)
    density.mixture2<-rowSums(densities.2)

    num<-(indicadoras*densities)
    den<-rowSums(num)
    prob_clusters<-num/den

    jj<-matrix(c(1:nrow(prob_clusters)),ncol=1)
    cluster.0<-apply(jj,1, function (jj) sample(k,1,replace=FALSE, prob=prob_clusters[jj,]))

    data.initial[,3]<-cluster.0


    # Updating the hypothesis

    if(hypothesis==1){

      joint.1<-function(l) dnorm(l[2], mean = mu.0[2], sd = sqrt(kappa), log = TRUE)+
        dgamma(1/l[4],shape=s,rate=s,log = TRUE)

      log.post.zeta.1<-log(pi.0)+sum(apply(parameters,2,joint.1))

      joint.0<-function(l) dnorm(l[2], mean = mu.0[2], sd = sqrt(kappa*v0), log = TRUE)+
        dgamma(1/l[4],shape=s*(1/v0),rate=s*(1/v0),log = TRUE)

      log.post.zeta.0<- log(1-pi.0)+sum(apply(parameters,2,joint.0))
    }

    if(hypothesis==0){

      joint.1<-function(l) dnorm(l[2], mean = mu.0[2], sd = sqrt(kappa*(1/v0)), log = TRUE)+
        dgamma(1/l[4],shape=s*v0,rate=s*v0,log = TRUE)

      log.post.zeta.1<-log(pi.0)+sum(apply(parameters,2,joint.1))

      joint.0<-function(l) dnorm(l[2], mean = mu.0[2], sd = sqrt(kappa), log = TRUE)+
        dgamma(1/l[4],shape=s,rate=s,log = TRUE)

      log.post.zeta.0<- log(1-pi.0)+sum(apply(parameters,2,joint.0))
    }

    vec.prob.zeta<-c(exp(log.post.zeta.1)/(exp(log.post.zeta.1)+exp(log.post.zeta.0)),exp(log.post.zeta.0)/(exp(log.post.zeta.1)+exp(log.post.zeta.0)))
    hypothesis<-sample(c(1,0),1,replace=F, prob=vec.prob.zeta)

    vec.zetas[a]<-hypothesis
    final.parameters[[a]]<-rbind(parameters,p.j)


    # Updating de Pi, Kappa and S

    if(hypothesis==1){

      pi.0<-rbeta(1,(1/2)+1,(3/2)-1)

      kappa<-1/rgamma(1, a2+(N/2),b2+(sum(parameters[2,]))^2/(2))

      log.post<-function(x){
        (a3-1)*log(x) - x*b3 + ncol(parameters)*x*log(x) - ncol(parameters)*(lgamma(x))+
          (x-1)*sum(log(1/parameters[4,])) - x*(sum(1/parameters[4,]))
      }

      stars.s2<-truncnorm::rtruncnorm(1, a=0.001, b=Inf, mean = s, sd = 0.3)

      log.r<-log.post(stars.s2)-log(truncnorm::dtruncnorm(stars.s2, a=0.001, b=Inf, mean = s, sd = 0.3))-
        log.post(s)+log(truncnorm::dtruncnorm(s, a=0.001, b=Inf, mean = stars.s2, sd = 0.3))

      if(log(runif(1))< min(0,log.r)){
        s<-stars.s2
      }

    }else{

      pi.0<-rbeta(1,1/2,3/2)

      kappa<-1/rgamma(1, a2+(N/2), b2+ (sum(parameters[2,]))^2/(2*v0))

      log.post<-function(x){
        (a3-1)*log(x/v0) - x*b3 + ncol(parameters)*(x/v0)*log(x/v0) - ncol(parameters)*(lgamma(x/v0))+
          ((x/v0)-1)*sum(log(1/parameters[4,])) - (x/v0)* (sum(1/parameters[4,]))
      }

      stars.s2<-truncnorm::rtruncnorm(1, a=0.001, b=Inf, mean = s, sd = 0.3)

      log.r<-log.post(stars.s2)-log(truncnorm::dtruncnorm(stars.s2, a=0.001, b=Inf, mean = s, sd = 0.3))-
        log.post(s)+log(truncnorm::dtruncnorm(s, a=0.001, b=Inf, mean = stars.s2, sd = 0.3))

      if(log(runif(1)) < min(0,log.r)){
        s<-stars.s2
      }

    }


    densities_1[a]<-list(density.mixture1)
    densities_2[a]<-list(density.mixture2)
  }


  burned.parameters<-final.parameters[-c(1:(n.mcm*0.2))]
  burned.hypothesis<-vec.zetas[-c(1:(n.mcm*0.2))]

  final.params<-burned.parameters[seq(from=1, to=length(burned.parameters), by=8)]
  final.hypo.vec<-burned.hypothesis[seq(from=1, to=length(burned.hypothesis), by=8)]

  post.probabilities<-sum(final.hypo.vec)/length(final.hypo.vec)

  final_list<-list(sampling.parameters=final.params,posterior.probability.H1=round(post.probabilities,5), data.init=sample.y)

  w<-matrix(1:n.mcm,ncol=1)
  den1<-apply(w,1,function(w) densities_1[[w]])
  den2<-apply(w,1,function(w) densities_2[[w]])

  burned.densities1<-den1[,-c(1:(n.mcm*0.2))]
  burned.densities2<-den2[,-c(1:(n.mcm*0.2))]

  m1<-seq(1, ncol(burned.densities1), by=8)
  final_densities1<-burned.densities1[,m1]
  final_densities2<-burned.densities2[,m1]


  # Posterior means ---------------------------------------------------------

  mean.posterior1<-rowSums(final_densities1)/ncol(final_densities1)
  mean.posterior2<-rowSums(final_densities2)/ncol(final_densities2)


  # Graph of marginal distributions -----------------------------------------

  matrix.q<-matrix(1:nrow(final_densities1),ncol=1)
  q1<-t(apply(matrix.q,1, function(j) quantile(final_densities1[j,], c(.025, .975))))
  q2<-t(apply(matrix.q,1, function(j) quantile(final_densities2[j,], c(.025, .975))))

  x3<-c(seq(min_value,max_value,length.out=200),rev(seq(min_value,max_value,length.out=200)))

  y31<-c(q1[,1],rev(q1[,2]))
  y32<-c(q2[,1],rev(q2[,2]))

  limit.y<-max(c(y31,y32))+0.1

  means_data_frame<-c(mean.posterior1,mean.posterior2)
  labels_means<-c(rep("E(g1(y)|data)", length(mean.posterior1)),rep("E(g2(y)|data)",length(mean.posterior2)))

  data.intervals <- data.frame(x3 = x3, y31 = y31, y32 = y32)
  data.means.estimations <- data.frame(means_data_frame,labels_means)

  `Grid value`<-round(c(rep(seq(min_value,max_value,length.out=200),2)),4)
  `Posterior density`<-round(means_data_frame,4)

  p<-ggplot2::ggplot()+geom_polygon(data=data.intervals, mapping=aes(x=x3, y=y31), fill = 'grey', colour = 'white') +
    geom_polygon(data=data.intervals, mapping=aes(x=x3, y=y32),fill = 'grey69', colour = 'white') +
    geom_line(data=data.means.estimations, mapping=aes(x=`Grid value`, y=`Posterior density`, color=labels_means))+
    labs(color="Posterior Means", x="y", y="g(y)") +scale_colour_manual(labels=c(expression(paste("E(", g[1](y),"|", "Data",")")),expression(paste("E(", g[2](y),"|", "Data",")"))), values=c("#330099","#993300")) + theme_bw()

  print(plotly::ggplotly(p, tooltip = c("x", "y")))


  if(post.probabilities<0.5){
    cat("\n There are no significant differences between marginal distributions G1 and G2.\n")
    print(paste("Posterior probability for H1:",round(post.probabilities,5)))
  }else{
    cat("\n There are significant differences between marginal distributions G1 and G2.\n")
    print(paste("Posterior probability for H1:",round(post.probabilities,5)))}

  final_list

}



#' @title Shift Plot Function
#'
#'
#' @description
#' The shift function quantifies the differences between the two study
#' populations, before and after the intervention or measurement, through
#' the difference between the quantiles of the two distributions as a
#' function of the quantiles of the first group.
#'
#'
#' @param results_BNP list of results obtained from the \code{BNP.test} function.
#'
#'
#' @return
#' It returns an interactive graph with the shift function estimate and the
#' credible interval at 95% confidence. By scrolling within the graph it is
#' possible to observe in detail the grid of 200 values ranging from the
#' minimum of the data to the maximum; it also shows the average of the
#' differences between the quantiles of both distributions.
#'
#'
#' @import plotly
#' @import ggplot2
#' @import dplyr
#'
#'
#' @note It is suggested to use the shift function if significant
#' differences are observed between the marginal distributions
#' (\code{posterior.probability.H1}), after performing the
#' hypothesis test \code{BNP.test}.
#'
#'
#' @examples
#'
#' \dontrun{
#' results <- BNP.test(x1, y1, n.mcm=10000)
#' plotshift.function(results)
#'}
#'
#' @export
plotshift.function <- function(results_BNP){

  data<-c(results_BNP$data.init[,1], results_BNP$data.init[,2])
  low.limit<-min(data)
  up.limit<-max(data)

  x.design<-matrix(c(1,1,0,1),ncol=2)

  # Cumulate marginal for Y1

  seq_marg<-seq(low.limit,up.limit,length.out=200)

  G1<-t(sapply(1:length(results_BNP[[1]]), function(i) rowSums(apply(results_BNP[[1]][[i]],2,function(x.data){
    x.data[8]*pnorm(q=seq_marg,mean=x.design[1,]%*%x.data[1:2],
                    sd=sqrt(x.data[3]+x.data[7]))})) ))

  shift.function<-function(R){

    # Cumulate marginal for Y2

    G2<-function(x){sum(apply(results_BNP[[1]][[R]],2,function(x.data,a){
      x.data[8]*pnorm(q=a,mean=x.design[2,]%*%x.data[1:2],
                      sd=sqrt((x.data[3]*x.data[4])+x.data[7]))},x))}

    # Finding the roots of the function

    inv.func = function (f, lower = -50000, upper = 50000) {
      function (y) uniroot((function (x) f(x) - y), lower = lower, upper = upper, extendInt = "yes")$root
    }

    q.inverse = inv.func(G2, -50000, 50000)

    matrix_k<-matrix(G1[R,],ncol=1)
    q.marginal.G2<-unlist(apply(matrix_k, 1, function(i) q.inverse(i)))
    shift.value<-q.marginal.G2-seq_marg
    return(shift.value)
  }

  R<-matrix(1:length(results_BNP$sampling.parameters),ncol=1)

  matrix_shift<-t(apply(R,1,shift.function))

  shift.function.means<-colMeans(matrix_shift)

  shift.function.quant<-t(apply(matrix_shift, 2, function(i) quantile(i, c(.025, .975))))

  x<-c(seq_marg,rev(seq_marg))

  y<-c(shift.function.quant[,1],rev(shift.function.quant[,2]))

  limit.y.low<-min(c(y))-0.1
  limit.y.up<-max(c(y))+0.1

  shift.data<-data.frame(x,y,shift.function.means)
  proof<-data.frame(seq_marg,shift.function.means)

  `Mean shift value`<-round(shift.function.means,4)
  `Grid value`<-round(seq_marg,4)

  hift<-ggplot2::ggplot(data=shift.data, aes(x=x,y=y, color=`Mean shift value`)) + geom_polygon( fill = 'grey', colour = 'white') +
    geom_line(data=proof, mapping=aes(x=`Grid value`, y=`Mean shift value`, colour="Shift estimation")) +
    theme(legend.position="top")+ scale_color_manual(name = "", values = c("Shift estimation" = "grey0")) +xlab('y')+ ylab('y2 - y1') + theme_bw()

  print(plotly::ggplotly(hift, tooltip = c("x", "y", "colour"))%>%
    plotly::layout(legend = list(orientation = "h", xanchor = "center", x = 0.5, y= 1.2)))
}



#' @title Contour Plot Function
#'
#'
#' @description
#'
#'The contour plot allows to visualize the joint distribution of the data,
#'as well as to identify the correlation, if any, between the two study populations.
#'
#'
#' @param results_BNP list of results obtained from the \code{BNP.test} function.
#'
#'
#' @return
#'
#' It returns an interactive graph showing the contours and the points representing
#' each pair of sample data. By scrolling within the graph it is possible to observe
#'  in detail the coordinates that delimit the contour and the density they represent.
#'
#'
#' @importFrom tidyr gather
#' @importFrom tibble rownames_to_column
#' @import dplyr
#' @import plotly
#' @import ggplot2
#'
#'
#' @examples
#'
#' \dontrun{
#' results <- BNP.test(x1, y1, n.mcm=10000)
#' contours.plot(results)
#'}
#'
#' @export
contours.plot <- function(results_BNP){

  data<-c(results_BNP$data.init[,1], results_BNP$data.init[,2])
  low.limit<-min(data)-1
  up.limit<-max(data)+1

  matrix_design<-matrix(c(1,1,0,1),ncol=2)

  x.coordinates<-seq(low.limit,up.limit,length.out=200)
  y.coordinates<-seq(low.limit,up.limit,length.out=200)

  density.function<-function(m, a, b){
    rowSums(apply(X=m,2,function(v,a,b){
      v[8]*mvtnorm::dmvnorm(x=cbind(a,b),mean=t(matrix_design%*%v[1:2]),
                            sigma=matrix(c(v[3]+v[7],rep((1-2*v[5])*(1-2*v[6])*v[7],2),
                                           (v[4]*v[3])+v[7]),ncol=2))
    },a,b))}

  z.coordinates<-matrix(rowMeans(sapply(X=results_BNP[[1]], FUN = density.function, a=expand.grid(x.coordinates,y.coordinates)$Var1,
                                          b=expand.grid(x.coordinates,y.coordinates)$Var2)),length(x.coordinates),length(y.coordinates))

  rownames(z.coordinates)<-seq(low.limit,up.limit,length.out=nrow(z.coordinates))

  colnames(z.coordinates)<-seq(low.limit,up.limit,length.out=ncol(z.coordinates))

  plot.contour<-as.data.frame(z.coordinates) %>%
    tibble::rownames_to_column() %>%
    tidyr::gather(y2, value, -rowname) %>%
    dplyr::mutate(y2 = as.numeric(y2),
           y1 = as.numeric(rowname)) %>%
    ggplot2::ggplot() +
    stat_contour(aes(x = y1, y = y2, z = value), color='grey0')+
    geom_point(data = data.frame(y1 = results_BNP$data.init[,1],
                                 y2 = results_BNP$data.init[,2]),
               mapping =aes(x = y1, y = y2), color='grey0', size=1)+labs(x='y1', y='y2') +
    theme_bw()

  print(plotly::ggplotly(plot.contour))
}
