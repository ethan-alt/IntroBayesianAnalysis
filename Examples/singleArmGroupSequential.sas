

%macro SA_SEQ(nSim=,pi0=,piActual=,nMax=,delta=,h=,alpha=,seed=1,out=_temp);
proc IML;
 call streaminit(&seed);
 nSim = &nSim;

 pi0     = &pi0;
 true_pi = &piActual;

 nMax    = &nMax;
 delta   = &delta;
 K0      = ceil(nMax/delta);
 h       = &h;
 alpha   = &alpha;


 alpha0 =    pi0 *K0*delta*h;
 gamma0 = (1-pi0)*K0*delta*h;

 x = t(do(0.001,0.999,0.001));
 y = pdf('beta',x,alpha0,gamma0);
 d = x||y;
 create plotdata from d[c={"pi" "pdf"}];
 	append from d;
 close plotdata;

 results = J(nSim,11,0);
 do sim = 1 to nSim;
	 y = 0;
	 n = 0;
	 do k = 1 to K0; 
		y = y + rand('binomial',true_pi,delta);
		n = n + delta;

		pp_Eff = sdf('beta',pi0,y+alpha0,(n-y)+gamma0);
		pp_Fut = cdf('beta',pi0,y+alpha0,(n-y)+gamma0);
		futility  = (pp_Eff<(alpha/2));
		efficacy  = (pp_Fut<(alpha/2));

		mean  = (y+alpha0)/(n+alpha0+gamma0);
		lower = quantile('beta',0.025,y+alpha0,(n-y)+gamma0);
		upper = quantile('beta',0.975,y+alpha0,(n-y)+gamma0);

		if (futility + efficacy) = 1 | k = K0 then do;
			results[sim,] = k||y||n||mean||lower||upper||futility||efficacy||(futility + efficacy)||pp_Eff||pp_Fut;
			k = K0 + 1;
		end;
	 end;
 end;

 create &out from results[c={"k" "y" "n" "mean" "lower" "upper" "futility" "efficacy" "error" "pp_Eff" "pp_Fut"}];
 	append from results;
 close &out;

quit;

 data &out;
  set &out;
   piActual = &piActual;
 run;

 data plotdata;
  set plotdata;
   low = 0;
 run;

 proc means data = &out maxdec=3 mean std median min max;
  var k n mean lower upper futility efficacy error pp_eff pp_fut;
  label k     = 'Number of Interim Analyses'
        n     = 'Sample Size'
        mean  = 'Sample Proportion'
        lower = 'Lower Limit of 95% CI'
		Upper = 'Upper Limit of 95% CI'
        futility = 'Futility'
        efficacy = 'Efficacy'
        error    = 'Futility or Efficacy'
        pp_eff   = 'Prob(Efficacy)'
        pp_fut   = 'Prob(Futility)';
 run;

proc sort data = &out out = futility;
	by k decending y;
	where futility = 1;
run;

data futility;
	set futility;
	by k;
	if first.k;
run;

proc sort data = &out out = efficacy;
	by k y;
	where efficacy = 1;
run;

data efficacy;
	set efficacy;
	by k;
	if first.k;
run;	

data combined;
 set futility(in=a)
 	 efficacy(in=b);
	 y=y/n;
	 if a then ylow = 0;
	 if b then yHigh  = 1;

	 if a then ya = y;
	 if b then yb = y;

run;

%mend;

title1 "Summary of Operating Characteristics Under Null Hypothesis";
%SA_SEQ(nSim=500000,pi0=0.30,piActual=0.30,
        nMax=90,delta=15,h=0.3,alpha=0.05,
        seed=1531,out=null);

title1 "Handicap Prior with (h=0.30)";
proc sgplot data = plotdata;
 	band x=pi lower=low upper=pdf;
	refline 0.3 / axis=x lineattrs=(pattern=2 thickness=1 color=lightred);
	xaxis label = 'Probability of Response';
	yaxis label = 'Density';
run;
quit;

title1 "Stopping Rules";
proc sgplot data = combined;
 	band x=n lower=yLow upper=y  / fillattrs=(color=verylightred)  nooutline legendlabel='Futility' name='A';
	band x=n lower=y upper=yHigh / fillattrs=(color=lightblue) nooutline legendlabel='Efficacy' name='B';
	scatter x=n y=ya / datalabel=mean markerattrs=(symbol=squareFilled color=red) ; format mean 6.2;
	scatter x=n y=yb / datalabel=mean markerattrs=(symbol=squareFilled color=blue) ; format mean 6.2;
	yaxis label = 'Proportion of Responders';
	xaxis label = 'Sample Size';
	keylegend "A" "B";
 run;
quit;


title1 "Summary of Operating Characteristics Under Alternative Hypothesis";
%SA_SEQ(nSim=500000,pi0=0.30,piActual=0.50,
        nMax=90,delta=15,h=0.3,alpha=0.05,
        seed=5423,out=alt);


