function [lambdaEstimationA, lambdaEstimationF, consensusA, consensusV, bertrand1A, bertrand1F, bertrand2A, bertrand2F, diLorenzoA, diLorenzoF, sahaiA] = errorInjection (A, F, consV, p, errormodel)
	NoisyA = makeAnoisy(A,p);
	failedLinks = A.-NoisyA;
	switch (errormodel)
		case 1
			lambdaEstimationA = NoisyA;
			lambdaEstimationF = F;
			consensusA = NoisyA;
			consensusV = consV;
			bertrand1A = NoisyA;
			bertrand1F = F;
	  		bertrand2A = NoisyA;
	  		bertrand2F = F;
	  		diLorenzoA = NoisyA;
	  		diLorenzoF = F;
	  		sahaiA = NoisyA;
		case 2
			lambdaEstimationA = NoisyA;
			lambdaEstimationF = F;
			consensusA = NoisyA;
			consensusV = consV;
			bertrand1A = NoisyA;
			bertrand1F = F;
	  		bertrand2A = NoisyA;
	  		bertrand2F = F;
	  		diLorenzoA = NoisyA;
	  		diLorenzoF = F;
	  		sahaiA = NoisyA;
		case 3
			#
		otherwise
			lambdaEstimationA = A;
			lambdaEstimationF = F;
			consensusA = A;
			consensusV = consV;
			bertrand1A = A;
			bertrand1F = F;
	  		bertrand2A = A;
	  		bertrand2F = F;
	  		diLorenzoA = A;
	  		diLorenzoF = F;
	  		sahaiA = A;
	endswitch
endfunction
