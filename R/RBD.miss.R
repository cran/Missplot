#' @export
#' @importFrom stats aov pf
#' @title Analysis of RBD when there is one missing observation
#' @param data A data-frame containing response in the first column , block number in the second column,
#' treatment number in the third column corresponding to the response value.
#' @param bl.miss Block number corresponding to the missing observation.
#' @param trt.miss Treatment number corresponding to the missing observation.
#' @return A data-frame containing x.hat , SSE x.hat , x_double.hat , SSE x_double.hat,
#' F statistics , p-value.
#' @details In design of experiments in RBD setup if there is one missing observation
#' present in the design , we can use the function RBD.miss to estimate the missing observation for testing
#' the differential effects for the treatments. Here, we estimate the missing obsevation by
#' minimizing the SSE of the design.
#' @examples
#' #Observation corresponding to the second block and third treatment is missing in the data
#' data=data.frame(res=rnorm(16,35,20),
#' block_no=rep(1:4,each=4),
#' trt_no=rep(1:4,times=4))
#' RBD.miss(data,2,3)
#' @author  Shantanu Nayek , Saheli Datta
#' @section Credits: Credits to Professor Surupa Chakraborty for building the theoritical concepts of Design of Experiment
#' and Professor Madhura Dasgupta for basic concepts for R.



RBD.miss=function(data,bl.miss,trt.miss)
{
  data
  y=data[,1]
  block=data[,2]
  treat=data[,3]
  b=length(unique(block))
  v=length(unique(treat))

  #Obtaining the missing observation
  miss.obs=y[block==bl.miss&treat==trt.miss];miss.obs

  #Computation of x.hat
  Bi=sum(y[block==bl.miss])-miss.obs;Bi
  Tj=sum(y[treat==trt.miss])-miss.obs;Tj
  G=sum(y)-miss.obs;G

  #Obtaining the missing observation
  x.hat=(b*Bi+v*Tj-G)/((b-1)*(v-1));x.hat


  y[block==bl.miss&treat==trt.miss]=x.hat
  y
  #Obtaining the SSE

  anova=summary(aov(y~as.factor(block)+as.factor(treat)))
  anova

  SSE.xhat=anova[[1]]$`Sum Sq`[3]
  df_SSE.xhat=anova[[1]]$Df[3]-1



  #Computation of missing observation under Ho
  xd.hat=Bi/(v-1)
  y[block==bl.miss&treat==trt.miss]=xd.hat
  y
  #Obtaining the sum of squares
  anovaHo=summary(aov(y~as.factor(block)))
  anovaHo

  SSEHo=anovaHo[[1]]$`Sum Sq`[2]
  df_SSEHo=anova[[1]]$Df[2]-1

  #Obtaining F statistics
  Fstat=((SSEHo-SSE.xhat)/(v-1))/(SSE.xhat/df_SSE.xhat);Fstat


  #p- value
  pv=pf(Fstat,v-1,df_SSE.xhat);pv

  output=data.frame(xhat=x.hat, SSExhat=SSE.xhat,
                    xdouble_hat=xd.hat ,SSEHo_xdhat=SSEHo,
                    Fstatistics=Fstat,pvalue=pv)
  output

}
