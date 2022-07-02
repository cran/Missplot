#' @export
#' @importFrom stats aov pf
#' @title Analysis of LSD when there is one missing observation
#' @param data A data-frame containing response in the first column ,row number in the second column,
#' column number in the third column , treatment number in the fourth column corresponding to the response value.
#' @param row.miss Row number corresponding to the missing observation.
#' @param col.miss Column number corresponding to the missing observation.
#' @return A data-frame containing x.hat , SSE x.hat , x_double.hat , SSE x_double.hat,
#' F statistics , p-value.
#' @details In design of experiments in LSD setup if there is one missing observation
#' present in the design , we can use the function LSD.miss to estimate the missing observation for testing
#' the differential effects for the treatments. Here, we estimate the missing obsevation by
#' minimizing the SSE of the design.
#' @examples
#' #Observation corresponding to the second row and third column is missing in the data
#' data=data.frame(res=rnorm(16,35,20),row_no=rep(1:4,each=4),col_no=rep(1:4,times=4),
#'     treat=c(1,2,3,4,2,3,4,1,3,4,1,2,4,1,2,3))
#' LSD.miss(data,2,3)
#' @author Saheli Datta , Shantanu Nayek
#' @section Credits: Credits to Professor Surupa Chakraborty for building the theoritical concepts of Design of Experiment
#' and Professor Madhura Dasgupta for basic concepts for R.
#' @section  Remark: Information on row number and column number corresponding to the missing observation
#' is to be known.



#Dataframe : Response, Row , Column , Treatment
LSD.miss=function(data,row.miss,col.miss)
{
  v=length(unique(data[,2]))
  #Obtaining the missing observation
  y=data[,1]
  row=data[,2]
  col=data[,3]
  trt=data[,4]

  miss.obs=y[row==row.miss&col==col.miss];miss.obs
  t.miss=trt[row==row.miss&col==col.miss];t.miss
  #Computation of x.hat
  Ri=sum(y[row==row.miss])-miss.obs;Ri
  Cj=sum(y[col==col.miss])-miss.obs;Cj
  Tk=sum(y[trt==as.factor(t.miss)])-miss.obs;Tk
  G=sum(y)-miss.obs;G

  x.hat=(v*(Ri+Cj+Tk)-2*G)/((v-1)*(v-2));x.hat
  y[row==row.miss&col==col.miss]=x.hat
  y
  #Obtaining the SSE

  anova=summary(aov(y~as.factor(row)+as.factor(col)+as.factor(trt)))
  anova

  SSE.xhat=anova[[1]]$`Sum Sq`[4]
  df_SSE=anova[[1]]$Df[4]-1

  #Computation of missing observation under Ho
  xd.hat=(v*(Ri+Cj)-G)/(v-1)^2;xd.hat
  y[row==row.miss&col==col.miss]=xd.hat
  y
  #Obtaining the sum of squares
  anovaHo=summary(aov(y~as.factor(row)+as.factor(col)))
  anovaHo

  SSEHo=anovaHo[[1]]$`Sum Sq`[3]
  df_SSEHo=anovaHo[[1]]$Df[3]-1

  #Obtaining F statistics
  Fstat=((SSEHo-SSE.xhat)/(v-1))/(SSE.xhat/df_SSE)
  #p- value
  pv=pf(Fstat,v-1,df_SSE);pv

  output=data.frame(xhat=x.hat, SSExhat=SSE.xhat,
                    xdouble_hat=xd.hat ,SSEHo_xdhat=SSEHo,
                    Fstatistics=Fstat,pvalue=pv)
  output


}

