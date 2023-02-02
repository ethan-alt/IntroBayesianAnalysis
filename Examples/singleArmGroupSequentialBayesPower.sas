

%macro SA_SEQ_BP(nSim=,pi0=,shape1=,shape2=,nMax=,delta=,h=,alpha=,seed=1,out=_temp);
proc IML;
 call streaminit(&seed);
 nSim = &nSim;

 pi0     = &pi0;
 nMax    = &nMax;
 delta   = &delta;
 K0      = ceil(nMax/delta);
 h       = &h;
 alpha   = &alpha;

 ** create hyperparameters for handicap prior;
 alpha0 =    pi0 *K0*delta*h;
 gamma0 = (1-pi0)*K0*delta*h;

 ** loop over simulation studies;
 results = J(nSim,6,0);
 do sim = 1 to nSim;
 	 true_pi = 0;
 	 do until(true_pi>pi0);
 	 	true_pi = rand('beta',&shape1,&shape2);
 	 end;
	 y       = 0;
	 n       = 0;
	 do k = 1 to K0; 
		y = y + rand('binomial',true_pi,delta);
		n = n + delta;

		pp_Eff    = sdf('beta',pi0,y+alpha0,(n-y)+gamma0);
		pp_Fut    = cdf('beta',pi0,y+alpha0,(n-y)+gamma0);
		futility  = (pp_Eff<(alpha/2));
		efficacy  = (pp_Fut<(alpha/2));
		mean      = (y+alpha0)/(n+alpha0+gamma0);


		if (futility + efficacy) = 1 | k = K0 then do;
			results[sim,] = k||y||n||mean||futility||efficacy;
			k = K0 + 1;
		end;
	 end;
 end;

 create &out from results[c={"k" "y" "n" "mean" "futility" "efficacy"}];
 	append from results;
 close &out;

quit;

 proc means data = &out maxdec=3 mean std median min max;
  var k n mean  futility efficacy;
  label k        = 'Number of Interim Analyses'
        n        = 'Sample Size'
        mean     = 'Sample Proportion'
        futility = 'Futility'
        efficacy = 'Efficacy';
 run;
%mend;

title1 "Summary of Operating Characteristics Under Alternative Hypothesis";
%SA_SEQ_BP(nSim=50000,pi0=0.30,shape1=50,shape2=50,
           nMax=90,delta=15,h=0.3,alpha=0.05,
           seed=1531);



