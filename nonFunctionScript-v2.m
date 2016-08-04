#! /usr/bin/octave
#Ubuntu's interpreter /usr/bin/octave
#OSX's interpreter /usr/local/bin/octave

addpath("~/Dropbox/papers/2015/aamas/");

	load ./matrices/A10s.mat; 
	#load ./matrices/A100s.mat; 
	#load ./matrices/A1000s.mat;

	#list parameters
	algorithmlist = [1]; #[1,2,3];
	swarmsizelist = [1]; #[1,2,3];
	packetdroplist = [0.01]; #,0.01, 0.05,0.1,0.2, 0.4,0.6,0.8];
	
	#scalar parameters
	numsamples = 10;
	maxiterations = 2001;
	errormodel = 1; #1: standard link failure with packet drop probability, 2: gentler error model, 3: ...

	fileid = 0;
	for algorithm=algorithmlist
		switch (algorithm)
	  	case 1 name = "bertrand";
	  	case 2 name = "dilorenzo";
	  	case 3 name = "sahai";
		endswitch
		for swarmsize=swarmsizelist
			switch (swarmsize)
		  	case 1 currentCell = A10s;
		  	case 2 currentCell = A100s;
		  	case 3 currentCell = A1000s;
			endswitch
			for packetdrop=packetdroplist
				switch (errormodel)
			  	case 1 errormodelname = "standard";
			  	case 2 errormodelname = "gentler";
			  	case 3 errormodelname = "tbd";
				endswitch
				fileid = fopen(strcat("./data/errorvssteps-",name,"-",int2str(swarmsize),"-",num2str(packetdrop),"-",errormodelname,".dat"), "w");
				fprintf(fileid, "x y errorx errory\n");
				avgtime = 0.0;
				printf(cstrcat("\n", name, " size=", int2str(swarmsize), " p=", num2str(packetdrop), " (", num2str(errormodel), ")\n"));
				#numsamples = length(currentCell);
				singledatapoint = zeros(numsamples,maxiterations);
				for sample=1:numsamples
					clk = tic;
					currentA = currentCell{sample};
					currentK = rows(currentA);
					[v,l] = eig(extractL(currentA));
					trueLambda2 = l(2,2);
					expectedOmega2 = preComputeOmega(currentA,maxiterations);
					###
					#algorithms variables
					bertrandPhase = 0; #0,1,2
					bertrandGcounter = 0;
					dilorenzoPhase = 0; #0,1,2
					lastF = ([0:currentK-1]'/currentK)./norm([0:currentK-1]'/currentK);
					consensusVector = lastF;
					newF = zeros(currentK,1);
					Us = zeros(maxiterations+2,currentK); Us(2,:) = unifrnd(0.0,1.0,[1,currentK]); Us(1,:) = Us(2,:);
					omegas = zeros(currentK,currentK);
					###
					for iteration=1:maxiterations
						###
						[lambdaEstimationA, lambdaEstimationF, consensusA, consensusV, bertrand1A, bertrand1F, bertrand2A, bertrand2F, diLorenzoA, diLorenzoF, sahaiA] = errorInjection (currentA, lastF, consensusVector, packetdrop, errormodel);
						bertrand1Degrees = sum(bertrand1A,2);
						bertrand2Degrees = sum(bertrand2A,2);
						diLorenzoDegrees = sum(diLorenzoA,2);
						###
						#compute/print error
						if algorithm == 1 || algorithm == 2
							estimatedLambda2s = (diag(extractL(lambdaEstimationA)) - ((lambdaEstimationA*lambdaEstimationF)./lambdaEstimationF))';
							differences = abs(estimatedLambda2s.-trueLambda2);
							err = (median(differences(isfinite(differences))))/trueLambda2;
							singledatapoint(sample,iteration) = err;
							printf(">sample %d (avg. time %f) iteration %d, err %f     \r", sample, avgtime, iteration, err);
						elseif algorithm == 3
							estimatedOmegas2s = omegas(:,2)';
							differences = abs(estimatedOmegas2s.-expectedOmega2);
							err = (median(differences(isfinite(differences))))/expectedOmega2;
							singledatapoint(sample,iteration) = err;
							printf(">sample %d (avg. time %f) iteration %d, err %f     \r", sample, avgtime, iteration, err);
						endif
						#run algorithm
					  	switch (algorithm)
						case 1
							if bertrandPhase == 0 #multiply by L
								#
								for i=1:currentK
									newF(i,1) = bertrand1Degrees(i,1)*bertrand1F(i,1)-bertrand1A(i,:)*bertrand1F;
								endfor
								lastF = newF;
								#
								bertrandPhase = 1; continue;
							elseif bertrandPhase == 1 #multiply by G
								#
								bertrandGcounter++;
								for i=1:currentK
									newF(i,1) = ((1-(bertrand2Degrees(i,1)/(3*currentK)))*bertrand2F(i,1))+(1/(3*currentK))*(bertrand2A(i,:)*bertrand2F);
								endfor
								lastF = newF;
								#
								if bertrandGcounter == 100
									bertrandGcounter = 0;
									bertrandPhase = 2; consensusVector = lastF; continue;
								endif
							elseif bertrandPhase == 2 #normalize with consenus average
								consensusVector = consensusIteration(consensusA,consensusV);
								if abs(max(consensusVector) - min(consensusVector)) < 0.01
									newF = newF./norm(newF);
									lastF = newF;
									bertrandPhase = 0; continue;
								endif
							endif
					  	case 2
							if dilorenzoPhase == 0 #share non-neighbors information with consenus average
								consensusVector = consensusIteration(consensusA,consensusV);
								if abs(max(consensusVector) - min(consensusVector)) < 0.01
									dilorenzoPhase = 1; continue;
								endif
							elseif dilorenzoPhase == 1 #multiply by B
								#
								#NOTE TO SELF: shouldn't I use the exact A to compute oppA as the previous consensus fase converged successfully?
								#
								oppA = abs(diLorenzoA-ones(currentK,currentK))-eye(currentK);
								for i=1:currentK
									one = (((currentK-1)/currentK) - diLorenzoDegrees(i,1)*(0.4*(2/currentK)))*diLorenzoF(i,1);
									two = ((currentK*(0.4*(2/currentK))-1)/currentK).*diLorenzoA(i,:)*diLorenzoF;
									three = (-1/currentK).*oppA(i,:)*diLorenzoF;
									newF(i,1) = one + two + three;
								endfor
								#
								lastF = newF;
								dilorenzoPhase = 2; consensusVector = newF; continue;
							elseif dilorenzoPhase == 2 #normalize with consenus average
								consensusVector = consensusIteration(consensusA,consensusV);
								if abs(max(consensusVector) - min(consensusVector)) < 0.01
									newF = newF./norm(newF);
									lastF = newF;
									dilorenzoPhase = 0; consensusVector = newF; continue;
								endif
							endif
					  	case 3
							#
							L = extractL(sahaiA);
							LNorm = zeros(currentK,currentK);
							for i=1:currentK
								if L(i,i) != 0
									LNorm(i,:) = L(i,:)./L(i,i);
								endif
							endfor
							for i=1:currentK
								Us(iteration+2,i) = 2*Us((iteration+2)-1,i) - Us((iteration+2)-2,i) - 0.8^2*(LNorm(i,:)*Us((iteration+2)-1,:)');
							endfor
							minimum = min([100,iteration+2]);
							Y = zeros(minimum,currentK);
							#Y = zeros(iteration+2,currentK);
							for i=1:currentK
								maximum = max([1,iteration-97]);
								Y(:,i) = fft(Us(maximum:iteration+2,i));
								#Y(:,i) = fft(Us(1:iteration+2,i));
							endfor
							omegas = zeros(currentK,currentK);
							for i=1:currentK
								len = rows(Y);
								[peak_list, indexes] = findpeaks(abs(Y(1:len/2,i)));
								omegas(i,1:length(indexes)) = 2*pi*(1/len).*(indexes.-1);
							endfor
							#
						endswitch
						###
					endfor
					avgtime = (avgtime*(sample-1) + toc(clk))/sample;
				endfor
				###
				#write files
				for iteration=1:maxiterations
					if mod(iteration, 50) == 1
						fprintf(fileid, "%f %f 0.0 %f\n", iteration, mean(singledatapoint(:,iteration)), std(singledatapoint(:,iteration)));
					endif
				endfor
				fclose(fileid);
				###
			endfor
		endfor
	endfor
