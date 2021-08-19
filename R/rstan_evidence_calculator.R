#'runs a the stan sampling of a model for some data and a vector of beta values
#'
#'@param model a rstan model
#'@param mydata a named list providing the data for the stan model
#'@param betavector a vector of values between 0 and 1 for beta
#'@param betalabel the label in the model that is used fot beta
runSampling<-function(model,mydata,betavector,betalabel="beta",...)
{
  run<-function(beta,model,data,betalabel,...)
  {
    print(list(beta=beta))
    cat("\n\n")
    data[betalabel]<-beta;
    rstan::sampling(model,data,diagnostic_file="diag",...);
  }
  f<-lapply(betavector,run,model=model,data=mydata,betalabel,...);
}

parameter_summary<-function(rstanrun, parameter)
{
  s=summary(rstanrun,parameter);
  mean=mean(s$c_summary[,"mean",]);
  sd_mean=sd(s$c_summary[,"mean",]);
  se_mean=sd_mean/(dim(s$c_summary)[3]-1)^0.5;
  sd=mean(s$c_summary[,"sd",]);
  sd_sd=sd(s$c_summary[,"sd",]);
  se_sd=sd_sd/(dim(s$c_summary)[3]-1)^0.5;
  var=mean(s$c_summary[,"sd",]^2);
  sd_var=sd(s$c_summary[,"sd",]^2);
  se_var=sd_var/(dim(s$c_summary)[3]-1)^0.5;
  c(mean=mean,se_mean=se_mean,sd=sd,se_sd=se_sd,var=var,se_var=se_var);
}

integrate_beta_with_derivative_mean_se<-function(beta,logL,se_logL,logL_var,se_logL_var)
{
  #    we assume that beta is sorted from 0 to 1.
  #
  # Obtains the Evidence from  runs of parallel models in a linear path between the prior and
  # the posterior likelihood.
  #  The formula of the probability function for each run
  #   logP(x,beta)=logPrior(x) + beta*logLikelihood(x), for 0=<beta=<1
  #
  #   this function takes as inputs parameters calculated after running such models. More
  #  precisely it accepts
  #   beta--> vector of the values of beta
  #   logL---> mean of the logLikelihood at each of beta values
  #   logL_se ---> error in the mean of logL
  #   dlogL---> mean of the variance of logL at each run
  #   dlogL_se---> error in the value of the variance of logL at each run
  #
  #   This function make use of the fact that the variance of logL is an estimator of its
  #   derivative with respect to beta. Then we have for each interval of beta the following
  #   values:

  #   At the start of the interval:   logL_0 and dlogL_0
  #   At the end of the interval:     logL_1 and dlogL_1
  #
  #   with the exception of the first and last value, each value works as the end of the
  #   previous interval and the begining of the net.
  #
  #   The trick to interpolate values of logL within the interval is the following:
  #   we use the values of logL_0 and logL_1 to calculate the average derivative within the
  #   interval and we divide the interval in two sections  of different second derivative of
  #   the logL. The first interval will go from the beginning to the deduced point where the
  #   dlogL reaches the mean dlogL value and the second interval from this point to the end.
  #   So we have two intervals  one where dlogL is above the mean and one below the mean, #
  #  therefore the integral of both intervals have to cancel each other. This fact allow us
  #  to find an equation to calculate the length of both sub-intervals
  #      delta_beta_0*(dlogL_0-dlogLmean)=delta_beta_1*(dlogLmean-dlogL1)
  #   and off course
  #           delta_beta_0+delta_beta_1=delta_beta
  #     then
  #       delta_beta_0= (dlogLmean-dlogL1)/(dlogL0-dlogL1)*delta_beta
  #
  #    with the value of both intervals we can calculate the second derivative of each
  #  interval and integrate this second order Taylor approximation to obtain the contribution
  #  to the Evidence of each subinterval.
  #
  #
  #   The calculation of the errors is straightforward using derivatives:
  #
  #      Ev=f(logL,dlogL)
  #
  #      se_Ev^2=(d_f/d_logL*se_logL)^2+(df/d_dlogL*se_logL)^2
  #
  #    We use the chain rule to propagate the errors from delta_beta_0/1 to logL_0/1 to logL
  #
  #    To estimate a supremum value for the approximation on the evolution of dlogL we
  #    consider that dlogL_0 remains constant during all the first interval and drops
  #    instantaneously to dlogL_1 for the all the second interval. This estimation is achieved
  #    by not taking into the second derivative of dlogL.
  #
  #    An infimun value would be obtained by assuming that the dlogL drops intantaneously to
  #   the mean value for all the interval. That value is obtained not using the values of
  #   dlogLs at all.
  #
  #    sometimes is not possible to calculate the value of logL for beta=0,
  #     so we add a correction for calculating the integral between 0 and the first value of
  #     beta, just by extrapolating the value of the second derivative.
  #
  #
  logL_se=se_logL;
  dlogL=logL_var;
  dlogL_se=se_logL_var;

  logL_0=logL[-length(logL)];
  logL_1=logL[-1];

  dlogL_0=dlogL[-length(dlogL)];
  dlogL_1=dlogL[-1];

  delta_beta<-diff(beta);

  delta_beta_0=(-dlogL_1 * delta_beta - logL_0+logL_1)/(dlogL_0-dlogL_1);

  d_delta_beta_0__d_logL_0= -1.0/(dlogL_0-dlogL_1);
  d_delta_beta_0__d_logL_1= 1.0/(dlogL_0-dlogL_1);
  d_delta_beta_0__d_dlogL_0= -(-dlogL_1 * delta_beta - logL_0+logL_1)/(dlogL_0-dlogL_1)^2;
  d_delta_beta_0__d_dlogL_1= (-dlogL_0 * delta_beta - logL_0+logL_1)/(dlogL_0-dlogL_1)^2;

  delta_beta_1=(dlogL_0 * delta_beta+logL_0-logL_1)/(dlogL_0-dlogL_1);
  d_delta_beta_1__d_logL_0= 1.0/(dlogL_0-dlogL_1);
  d_delta_beta_1__d_logL_1= -1.0/(dlogL_0-dlogL_1);
  d_delta_beta_1__d_dlogL_0= -(dlogL_1 * delta_beta + logL_0-logL_1)/(dlogL_0-dlogL_1)^2;
  d_delta_beta_1__d_dlogL_1=(dlogL_0 * delta_beta+logL_0-logL_1)/(dlogL_0-dlogL_1)^2;

  Ev_0=delta_beta_0*logL_0+0.5*delta_beta_0^2*(2/3*dlogL_0+1/3*(logL_1-logL_0)/delta_beta);
  d_Ev_0__d_delta_beta_0=logL_0+delta_beta_0*(2/3*dlogL_0+1/3*(logL_1-logL_0)/delta_beta);

  d_Ev_0__d_logL_0=delta_beta_0+
    -1/6*delta_beta_0^2/delta_beta+
    d_Ev_0__d_delta_beta_0*d_delta_beta_0__d_logL_0;

  d_Ev_0__d_logL_1=1/6*delta_beta_0^2/delta_beta+
    d_Ev_0__d_delta_beta_0*d_delta_beta_0__d_logL_1;

  d_Ev_0__d_dlogL_0=1/3*delta_beta_0^2+
    d_Ev_0__d_delta_beta_0*d_delta_beta_0__d_dlogL_0;


  d_Ev_0__d_dlogL_1=d_Ev_0__d_delta_beta_0*d_delta_beta_0__d_dlogL_1;



  Ev_1=delta_beta_1*logL_1-0.5*delta_beta_1^2*(2/3*dlogL_1+1/3*(logL_1-logL_0)/delta_beta);

  d_Ev_1__d_delta_beta_1=logL_1-delta_beta_1*(2/3*dlogL_1+1/3*(logL_1-logL_0)/delta_beta);

  d_Ev_1__d_logL_0=1/6*delta_beta_1^2/delta_beta+
    d_Ev_1__d_delta_beta_1*d_delta_beta_1__d_logL_0;
  d_Ev_1__d_logL_1=delta_beta_1-
    1.0/6.0*(delta_beta_1^2) /delta_beta+
    d_Ev_1__d_delta_beta_1*d_delta_beta_1__d_logL_1;
  d_Ev_1__d_dlogL_0=d_Ev_1__d_delta_beta_1*d_delta_beta_1__d_dlogL_0;
  d_Ev_1__d_dlogL_1=-1/3*delta_beta_1^2+
    d_Ev_1__d_delta_beta_1*d_delta_beta_1__d_dlogL_1;

  d_Ev__d_logL_0=d_Ev_0__d_logL_0+d_Ev_1__d_logL_0;
  d_Ev__d_logL_1=d_Ev_0__d_logL_1+d_Ev_1__d_logL_1;
  d_Ev__d_dlogL_0=d_Ev_0__d_dlogL_0+d_Ev_1__d_dlogL_0;
  d_Ev__d_dlogL_1=d_Ev_0__d_dlogL_1+d_Ev_1__d_dlogL_1;

  d_Ev__d_logL=c(d_Ev__d_logL_0,0)+c(0,d_Ev__d_logL_1);
  d_Ev__d_dlogL=c(d_Ev__d_dlogL_0,0)+c(0,d_Ev__d_dlogL_1);


  dlogLmean_subzero=(logL_1[1]-logL_0[1])/delta_beta[1];
  d2logL_subzero=(dlogL_0[1]-dlogLmean_subzero)/(dlogLmean_subzero-dlogL_1[1])*
    (dlogL_0[1]-dlogL_1[1])/delta_beta[1];

  d_dlogLmean_subzero__d_logL_1=1.0/delta_beta_1[1];
  d_dlogLmean_subzero__d_logL_0=-1.0/delta_beta_1[1];


  d_d2logL_subzero__d_dlogLmean_subzero=
    -(dlogL_0[1]-dlogL_1[1])^2/(dlogLmean_subzero-dlogL_1[1])^2/delta_beta[1];

  d_d2logL_subzero__d_dlogL_0=
    (1)/(dlogLmean_subzero-dlogL_1[1])*                (dlogL_0[1]-dlogL_1[1])/delta_beta[1]+
    (dlogL_0[1]-dlogLmean_subzero)/(dlogLmean_subzero-dlogL_1[1])/delta_beta[1];

  d_d2logL_subzero__d_dlogL_1=
    (dlogL_0[1]-dlogLmean_subzero)/(dlogLmean_subzero-dlogL_1[1])^2*
    (dlogL_0[1]-dlogL_1[1])/delta_beta[1]-
    (dlogL_0[1]-dlogLmean_subzero)/(dlogLmean_subzero-dlogL_1[1])/delta_beta[1];




  Ev_subzero=beta[1]*logL_0[1]-0.5*beta[1]^2*dlogL_0[1]-1/6*beta[1]^3*d2logL_subzero;

  d_Ev_subzero__d_logL_0=beta[1]-1/6*beta[1]^3*
    d_d2logL_subzero__d_dlogLmean_subzero*d_dlogLmean_subzero__d_logL_0;

  d_Ev_subzero__d_logL_1=-1/6*beta[1]^3*
    d_d2logL_subzero__d_dlogLmean_subzero*d_dlogLmean_subzero__d_logL_1;

  d_Ev_subzero__d_dlogL_0=-0.5*beta[1]^2-1/6*beta[1]^3*d_d2logL_subzero__d_dlogL_0;

  d_Ev_subzero__d_dlogL_1=-1/6*beta[1]^3*d_d2logL_subzero__d_dlogL_1;

  d_Ev__d_logL[1]<-d_Ev__d_logL[1]+d_Ev_subzero__d_logL_0;
  d_Ev__d_logL[2]<-d_Ev__d_logL[2]+d_Ev_subzero__d_logL_1;

  d_Ev__d_dlogL[1]<-d_Ev__d_dlogL[1]+d_Ev_subzero__d_dlogL_0;
  d_Ev__d_dlogL[2]<-d_Ev__d_dlogL[2]+d_Ev_subzero__d_dlogL_1;

  var_Ev_logL=(d_Ev__d_logL*logL_se)^2;
  var_Ev_dlogL=(d_Ev__d_dlogL*dlogL_se)^2;

  Ev_by_interval=c(Ev_subzero,Ev_0+Ev_1);
  Ev=sum(Ev_by_interval);
  se_Ev=sum(var_Ev_logL+var_Ev_dlogL)^0.5;

  Ev_0_supremum=delta_beta_0*logL_0+0.5*delta_beta_0^2*dlogL_0;
  Ev_1_supremum=delta_beta_1*logL_1-0.5*delta_beta_1^2*dlogL_1;

  Ev_subzero_supremum=beta[1]*logL_0[1]-0.5*beta[1]^2*dlogL_0[1];

  Ev_supremum_by_interval=c(Ev_subzero_supremum,Ev_0_supremum+Ev_1_supremum);

  Ev_error_above_by_interval=Ev_supremum_by_interval-Ev_by_interval;


  Ev_supremum=sum(Ev_supremum_by_interval);


  Ev_subzero_infimum=beta[1]*logL_0[1]-beta[1]^2*dlogL_0[1]-1/2*beta[1]^3*d2logL_subzero;



  Ev_infimum_by_interval=c(Ev_subzero_infimum,delta_beta*(logL_0+logL_1)/2);


  Ev_error_below_by_interval=Ev_infimum_by_interval-Ev_by_interval;


  Ev_infimum=sum(Ev_infimum_by_interval);



  res<-list(Ev=Ev,
            se_Ev=se_Ev,
            Ev_supremum=Ev_supremum,
            Ev_infimum=Ev_infimum,
            se_Ev_beta=rbind(beta=beta,
                             logL_se=var_Ev_logL^0.5,
                             dlogL_se=var_Ev_dlogL^0.5,
                             supremum=Ev_error_above_by_interval,
                             infimum=Ev_error_below_by_interval)
  );
  res
}
Evidencerun<-function(rstanrunvector,betavector,logLiklabel="loglikelihood")
{
  logL<-lapply(rstanrunvector,parameter_summary,parameter=logLiklabel);
  list(logL=logL,I=integrate_beta_with_derivative_mean_se(beta=betavector,
                                                          logL=map_dbl(logL, function(x) x["mean"]),
                                                          se_logL =map_dbl(logL, function(x) x["se_mean"]),
                                                          logL_var = map_dbl(logL, function(x) x["var"]),
                                                          se_logL_var=map_dbl(logL, function(x) x["se_var"])));

}




