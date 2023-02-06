%let seed  = 100;
%let nSim  = 100000;

data design_sim;
	call streaminit(&seed);
	n1 = 20 + 76/2;
	n0 = 4  + 76/2;

	do delta = 0.00 to 0.15 by 0.03;
	do a0    = 0.00 to 0.50 by 0.005;

		/** normal approximation for power prior **/
		shape01 = a0*((93  + 125))               + 0.5;
		shape02 = a0*((275 + 287) - (93  + 125)) + 0.5;

		shape11 = a0*((118 + 167))               + 0.5;
		shape12 = a0*((273 + 290) - (118 + 167)) + 0.5;

		pp1     = shape11/(shape11+shape12);
		pp0     = shape01/(shape01+shape02);

		ppMean  = pp1 - pp0;
		ppVar   = pp1*(1-pp1)/(shape11+shape12) + pp0*(1-pp0)/(shape01+shape02);

		tailArea = cdf('normal',0,ppMean,sqrt(ppVar));

        rej      = 0;
		meanDiff = 0;
		sqError  = 0;
		aveSD    = 0;
		badRej   = 0;
		do sim = 1 to &nsim;
			pi0   = rand('beta',(93  + 125),(275 + 287) - (93  + 125));
			pi1   = max(min(1,pi0+delta),0);

			y0    = rand('binomial',pi0,n0);	
			y1    = rand('binomial',pi1,n1);

            /** normal aproximation from likelihood **/
			np1    = y1/n1;
			np0    = y0/n0;

			npMean = np1 - np0;
			npVar  = np1*(1-np1)/n1 + np0*(1-np0)/n0;

			postMean = (npMean/npVar + ppMean/ppVar) /
			           (1/npVar      +      1/ppVar);
			postSD   = sqrt(1 / (1/npVar  +  1/ppVar));

			meanDiff = meanDiff + postMean;
			rej      = rej      + (sdf('normal',0,postMean,postSD)>0.975);
			sqError  = sqError  + (postMean-delta)**2;
			aveSD    = aveSD    + postSD;
			badRej   = badRej   + (sdf('normal',0,postMean,postSD)>0.975)*(y0>=y1);
		end;
		aveSD    = aveSD        / &nSim;
		meanDiff = meanDiff     / &nSim;
		sqError  = sqrt(sqError / &nSim);
		rej      = rej          / &nSim;
		badRej   = badRej       / &nSim;

		if round(a0,0.001) = 0.50 then deltac = put(delta,6.2);
		else deltac=' ';
		output;
	end;
	end;
	drop pi: pp: y1 y0 n: shape:;
run;

ods graphics / height=4in width=6in;
title "Mean Squared Error of Posterior Mean Estimator";
proc sgplot data = Design_Sim;
 loess x=a0 y=sqError / group=delta datalabel=deltac nomarkers;
 yaxis label = 'sqrt(MSE)';
 refline 0.098 / axis=y lineattrs=(pattern=2);
 refline 0.40   / axis=x lineattrs=(pattern=2);
run;


title "Null Hypothesis Rejection Rate";
proc sgplot data = Design_Sim;
 loess x=a0 y=rej / group=delta nomarkers;
 refline  0.34 0.80 / axis=y  axis=y lineattrs=(pattern=2);
 refline  0.26 / axis=x lineattrs=(pattern=2);
 yaxis label = 'Null Hypothesis Rejection Rate';
run;


