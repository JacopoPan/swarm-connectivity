#! /usr/local/bin/octave
#Ubuntu's interpreter /usr/bin/octave
#OSX's interpreter /usr/local/bin/octave

addpath("~/Dropbox/papers/2015/aamas/");

load ./matrices/A10s.mat; 
#load ./matrices/A100s.mat; 
#load ./matrices/A1000s.mat;

stepsfiles = zeros(3,3);
errorfiles = zeros(3,3);
for algorithm=[2]#[1,2,3]
	switch (algorithm)
  	case 1 name = "bertrand";
  	case 2 name = "dilorenzo";
  	case 3 name = "sahai";
	endswitch
	for swarmsize=[1]#[1,2,3]
		switch (swarmsize)
	  	case 1 currentCell = A10s;
	  	case 2 currentCell = A100s;
	  	case 3 currentCell = A1000s;
		endswitch
		stepsfiles(algorithm,swarmsize) = fopen(strcat("./data/stepsvsp-",name,"-",int2str(swarmsize),".dat"), "w");
		fprintf(stepsfiles(algorithm,swarmsize), "x y errorx errory\n");
		errorfiles(algorithm,swarmsize) = fopen(strcat("./data/errorvssteps-",name,"-",int2str(swarmsize),".dat"), "w");
		fprintf(errorfiles(algorithm,swarmsize), "x y errorx errory\n");
		for packetdrop=[0.0]#[0.0,0.4,0.8]#[0.0,0.2,0.4,0.6,0.8,0.95]
			maxiterations = 601;
			avgtime = 0.0;
			printf(cstrcat("\n", name, " size=", int2str(swarmsize), " p=", num2str(packetdrop), "\n"));
			#numsamples = length(currentCell);
			numsamples = 20;
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
				###
				printf(cstrcat(">sample=", int2str(sample)," avg time=", num2str(avgtime),"          \r"));
				for iteration=1:maxiterations
					currentNoisyA = makeAnoisy(currentA,packetdrop);
					#currentNoisyA = currentA;
					currentDegrees = sum(currentNoisyA,2);
					###
					#run algorithm
				  	switch (algorithm)
					case 1
						if bertrandPhase == 0 #multiply by L
							#
							for i=1:currentK
								newF(i,1) = currentDegrees(i,1)*lastF(i,1)-currentNoisyA(i,:)*lastF;
							endfor
							lastF = newF;
							#
							estimatedLambda2s = (diag(extractL(currentNoisyA)) - ((currentNoisyA*newF)./newF))';
							differences = abs(estimatedLambda2s.-trueLambda2);
							err = (median(differences(isfinite(differences))))/trueLambda2;
							singledatapoint(sample,iteration) = err;
							printf("iteration %d, err %f     \r", iteration, err); #input("aaa");
							bertrandPhase = 1; continue;
						elseif bertrandPhase == 1 #multiply by G
							#
							bertrandGcounter++;
							for i=1:currentK
								newF(i,1) = ((1-(currentDegrees(i,1)/(3*currentK)))*lastF(i,1))+(1/(3*currentK))*(currentNoisyA(i,:)*lastF);
							endfor
							lastF = newF;
							#
							estimatedLambda2s = (diag(extractL(currentNoisyA)) - ((currentNoisyA*newF)./newF))';
							differences = abs(estimatedLambda2s.-trueLambda2);
							err = (median(differences(isfinite(differences))))/trueLambda2;
							singledatapoint(sample,iteration) = err; 
							printf("iteration %d, err %f     \r", iteration, err); #input("aaa");
							if bertrandGcounter == 100
								bertrandGcounter = 0;
								bertrandPhase = 2; consensusVector = lastF; continue;
							endif
						elseif bertrandPhase == 2 #normalize with consenus average
							consensusVector = consensusIteration(currentNoisyA,consensusVector);
							if abs(max(consensusVector) - min(consensusVector)) < 0.01
								newF = newF./norm(newF);
								lastF = newF;
								estimatedLambda2s = (diag(extractL(currentNoisyA)) - ((currentNoisyA*lastF)./lastF))';
								differences = abs(estimatedLambda2s.-trueLambda2);
								err = (median(differences(isfinite(differences))))/trueLambda2;
								singledatapoint(sample,iteration) = err;
								printf("iteration %d, err %f     \r", iteration, err); #input("aaa");
								bertrandPhase = 0; continue;
							endif
							estimatedLambda2s = (diag(extractL(currentNoisyA)) - ((currentNoisyA*lastF)./lastF))';
							differences = abs(estimatedLambda2s.-trueLambda2);
							err = (median(differences(isfinite(differences))))/trueLambda2;
							singledatapoint(sample,iteration) = err;
							printf("iteration %d, err %f     \r", iteration, err); #input("aaa");
						endif
				  	case 2
						if dilorenzoPhase == 0 #share non-neighbors information with consenus average
							estimatedLambda2s = (diag(extractL(currentNoisyA)) - ((currentNoisyA*lastF)./lastF))';
							differences = abs(estimatedLambda2s.-trueLambda2);
							err = (median(differences(isfinite(differences))))/trueLambda2;
							singledatapoint(sample,iteration) = err;
							printf("iteration %d, err %f     \r", iteration, err); #input("aaa");
							consensusVector = consensusIteration(currentNoisyA,consensusVector);
							if abs(max(consensusVector)- min(consensusVector)) < 0.01
								dilorenzoPhase = 1; continue;
							endif
						elseif dilorenzoPhase == 1 #multiply by B
							#
							oppA = abs(currentNoisyA-ones(currentK,currentK))-eye(currentK);
							for i=1:currentK
								one = (((currentK-1)/currentK) - currentDegrees(i,1)*(0.4*(2/currentK)))*lastF(i,1);
								two = ((currentK*(0.4*(2/currentK))-1)/currentK).*currentNoisyA(i,:)*lastF;
								three = (-1/currentK).*oppA(i,:)*lastF;
								newF(i,1) = one + two + three;
							endfor
							#
							lastF = newF;
							estimatedLambda2s = (diag(extractL(currentNoisyA)) - ((currentNoisyA*lastF)./lastF))';
							differences = abs(estimatedLambda2s.-trueLambda2);
							err = (median(differences(isfinite(differences))))/trueLambda2;
							singledatapoint(sample,iteration) = err;
							printf("iteration %d, err %f     \r", iteration, err); #input("aaa");
							dilorenzoPhase = 2; consensusVector = newF; continue;
						elseif dilorenzoPhase == 2 #normalize with consenus average
							estimatedLambda2s = (diag(extractL(currentNoisyA)) - ((currentNoisyA*lastF)./lastF))';
							differences = abs(estimatedLambda2s.-trueLambda2);
							err = (median(differences(isfinite(differences))))/trueLambda2;
							singledatapoint(sample,iteration) = err;
							printf("iteration %d, err %f     \r", iteration, err); #input("aaa");
							consensusVector = consensusIteration(currentNoisyA,consensusVector);
							if abs(max(consensusVector)- min(consensusVector)) < 0.01
								newF = newF./norm(newF);
								lastF = newF;
								dilorenzoPhase = 0; consensusVector = newF; continue;
							endif
						endif
				  	case 3
						#
						L = extractL(currentNoisyA);
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
						omegas(:,2)';
						estimatedOmegas2s = omegas(:,2)';
						differences = abs(estimatedOmegas2s.-expectedOmega2);
						err = (median(differences(isfinite(differences))))/expectedOmega2;
						singledatapoint(sample,iteration) = err;
						printf("iteration %d, err %f     \r", iteration, err); #input("aaa");
						#
					endswitch
					###
				endfor
				avgtime = (avgtime*(sample-1) + toc(clk))/sample;
			endfor
			###
			#write files
			iterationstoconvergence = zeros(numsamples,1);
			for i=1:numsamples
				for j=1:maxiterations
					if (singledatapoint(i,j)) < 0.05 || j == maxiterations
						iterationstoconvergence(i,1) = j;
						break;
					endif
				endfor
			endfor
			fprintf(stepsfiles(algorithm,swarmsize), "%f %f 0.0 %f\n", packetdrop, median(iterationstoconvergence), max(iterationstoconvergence)-median(iterationstoconvergence));
			for iteration=1:maxiterations
				if mod(iteration, 100) == 1 && packetdrop == 0.0
					fprintf(errorfiles(algorithm,swarmsize), "%f %f 0.0 %f\n", iteration, mean(singledatapoint(:,iteration)), std(singledatapoint(:,iteration)));
				endif
			endfor
			###
		endfor
		fclose(stepsfiles(algorithm,swarmsize));
		fclose(errorfiles(algorithm,swarmsize));
	endfor
endfor