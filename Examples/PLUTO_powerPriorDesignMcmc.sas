%let seed  =   100;
%let nSim  =  2500;
%let nmc   = 10000;

data design_sim;
	n1 = 20 + 76/2;
	n0 = 4  + 76/2;

    /** these can be changed to a range **/
	do delta = 0.00 to 0.12 by 0.030;
	do a0    = 0.00 to 0.26 by 0.065;
	 output;
	end;
	end;
run;

data design_sim2;
 set design_sim;
  call streaminit(&seed.);
  
  do simulation = 1 to &nSim;
	  pi0   = rand('beta',(93  + 125),(275 + 287) - (93  + 125));
	  pi1   = max(min(1,pi0+delta),0);  
  
  	  source = 1;
	  weight = a0;
	  trt    = 0;
	  y      = (93  + 125);
	  n      = (275 + 287);
	   output;
	   
	  weight = a0;
	  trt    = 1;
	  y      = (118 + 167);
	  n      = (273 + 290);
	   output;   
	
	  source = 2;
	  weight = 1;
	  trt    = 0;
	  n      = n0;
	  y      = rand('binomial',pi0,n0);
	  output;

	  weight = 1;
	  trt    = 1;
	  n      = n1;
	  y      = rand('binomial',pi1,n1);
	  output;
	  
	end;
run;
proc sort; by delta a0 simulation; run;


option nonotes;
ods _all_ close;
ods output PostSumInt = PSI;
proc mcmc data = design_sim2 monitor=(diff pRej) nmc=&nmc.;
	by delta a0 simulation;
 
	 array pi[2] pi0c pi0h;
	 parms pi0c 0.35 pi0h 0.35 diff 0 / slice;
	  
	  prior pi0c ~ beta(0.5,0.5);
	  prior pi0h ~ beta(0.5,0.5);
	  
	  logPriorDiff = 0;
	  if not(0<pi0c+diff<1) then logPriorDiff = .;
	  if not(0<pi0h+diff<1) then logPriorDiff = .;
	  prior diff ~ general(logPriorDiff); 
	     
	  mu = pi0c*(source=2) + pi0h*(source=1) + diff*(trt=1);
	  
	  logLike = weight*logpdf('binomial',y,mu,n);
	  
	  beginnodata;
	   pRej  = (diff>0);
	  endnodata;
      
  model general(logLike);
quit;

options notes;
ods html newfile=none ;

data Results;
 set PSI;
  stat  = 'MSE  ';
  value = (mean-delta)**2;
  if parameter = 'diff' then output;
 
  stat  = 'Power';
  value = (mean>0.975);
  if parameter = 'pRej' then output;
run;

proc means data = Results nway noprint;
 class delta a0 stat;
 var value;
 output out = summary mean=value;
run;

data summary2; 
 set  summary;
 if stat = 'MSE' then value = sqrt(Value);
run;

ods pdf file="/home/map33232/short-course/programs/SAS/ex10-mcmc.pdf";

ods graphics / height=4in width=6in;
title "Mean Squared Error of Posterior Mean Estimator";
proc sgplot data = summary2;
 where stat = 'MSE';
 series x=a0 y=value / group=delta;
 yaxis label = 'sqrt(MSE)';
 refline 0.098 / axis=y lineattrs=(pattern=2);
 refline 0.40   / axis=x lineattrs=(pattern=2);
run;


title "Null Hypothesis Rejection Rate";
proc sgplot data = summary2;
 where stat = 'Power';
 series x=a0 y=value / group=delta;
 refline  0.34 0.80 / axis=y  axis=y lineattrs=(pattern=2);
 refline  0.26 / axis=x lineattrs=(pattern=2);
 yaxis label = 'Null Hypothesis Rejection Rate';
run;

ods pdf close;
