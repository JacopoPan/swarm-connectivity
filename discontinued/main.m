#! /usr/local/bin/octave
#Ubuntu's interpreter /usr/bin/octave
#OSX's interpreter /usr/local/bin/octave

addpath("~/Dropbox/papers/2015/aamas/");
load data.mat;
Z

N = 10;
T = 3;
l = 10.0;
r = 4.0;

positions = unifrnd(0,l,[N,2]);
A = zeros(N, N);
L = zeros(N,N);

function [newPositions] = moveRobots (oldPositions)
	N = rows(oldPositions);
	moves = unifrnd(-0.25,0.25,[N,2]);
	newPositions = oldPositions + moves;
	for i=1:N
		if newPositions(i,1) > 10.0 newPositions(i,1) = 10.0; endif
		if newPositions(i,1) < 0.0 newPositions(i,1) = 0.01; endif
		if newPositions(i,2) > 10.0 newPositions(i,2) = 10.0; endif
		if newPositions(i,2) < 0.0 newPositions(i,2) = 0.01; endif
	endfor
endfunction

function printTopology (positions)
	N = rows(positions);	
	printf("Topology\n");
	for i=1:10
		for j=1:30
			flag = 0;
			for k=1:N
				if (ceil(positions(k,1)*3) == j) && (ceil(positions(k,2)) == 11-i)
				printf("%d", k); flag = 1;
				endif
			endfor
			if !flag
			printf(" ");
			endif
		endfor
		printf("\n");
	endfor
	printf("\n");
endfunction

function newA = updateA (positions, r)
	N = rows(positions);
	newA = zeros(N,N);
	for i=1:(N-1)
		for j=(i+1):N
			if (norm(positions(i,:).-positions(j,:)) <= r)
			newA(i,j) = 1;
			newA(j,i) = 1;
			endif
		endfor
	endfor
endfunction

function L = extractL (A)
	N = rows(A);
	L = zeros(N,N);
	for i=1:N
		for j=i:N
			if i == j
			L(i,j) = sum(A(i,:),2);
			elseif A(i,j) == 1
			L(i,j) = -1;
			L(j,i) = -1;
			endif
		endfor
	endfor
endfunction

fid = zeros(N,3);
for i=1:N
	filename = strcat("plot", int2str(i), ".tex");
	id = fopen(filename, "w");
	fputs(id, strcat("\\begin{tikzpicture} \n \\begin{axis}[ \n height=3.5cm, \n width=18cm, ymin = -0.5, \n ymax = 0.5, \n grid=major, \n ] \n \\input{plot", int2str(i), "real.tex}\n \\input{plot", int2str(i), "cent.tex}\n \\input{plot", int2str(i), "dist.tex}\n \\end{axis}\n \\end{tikzpicture}"));
	fclose(id);

	for j=1:3
		chart = "";
		if j == 1 chart = "real";
		elseif j == 2 chart = "cent";
		else chart = "dist";
		endif
		filename = strcat("plot", int2str(i), chart, ".tex");
		fid(i,j) = fopen(filename, "w");
		fputs(fid(i,j), "\\addplot coordinates { \n");
	endfor
endfor

#parameters/initializations
centralizedPIxs = unifrnd(0.0,1.0,[N,1]);
centralizedDiLorenzo = unifrnd(0.0,1.0,[N,1]);
distrXs = unifrnd(0.0,1.0,[N,1]);
distrBs = unifrnd(0.0,1.0,[N,1]);
distrCs = unifrnd(0.0,1.0,[N,1]);
Tmax = 98;
c = 0.8;
waveUs = zeros(Tmax+2, N); waveUs(2,:) = unifrnd(0.0,1.0,[1,N]); waveUs(1,:) = waveUs(2,:);

iter = 0;
do
	iter = iter+1
	bakPositions = positions;
	#positions = moveRobots(positions); #whether we assume the robots are moving or not
	#A = updateA(positions,r); 
	#L = extractL(A); [v,D] = eig(L); v; lambdas = diag(D)';
	#if lambdas(1,2) < 0.001 && iter == 1
	#	do
	#	positions = unifrnd(0,l,[N,2]);
	#	A = updateA(positions,r); L = extractL(A); [v,D] = eig(L); v; lambdas = diag(D)';
	#	until lambdas(1,2) > 0.001
	#elseif lambdas(1,2) < 0.001
	#	positions = bakPositions; #refuse moves if they led to a disconnect
	#endif
	#printTopology(positions)
	#A = updateA(positions,r);
	A = [0,   1,   0,   0,   0,   1,   0,   0,   0,   0; 1,   0,   1,   1,   1,   1,   1,   1,   1,   1; 0,   1,   0,   0,   0,   0,   0,   1,   1,   1; 0,   1,   0,   0,   1,   0,   1,   0,   0,   0; 0,   1,   0,   1,   0,   0,   0,   0,   0,   0; 1,   1,   0,   0,   0,   0,   1,   0,   0,   0; 0,   1,   0,   1,   0,   1,   0,   0,   0,   0; 0,   1,   1,   0,   0,   0,   0,   0,   1,   1; 0,   1,   1,   0,   0,   0,   0,   1,   0,   0; 0,   1,   1,   0,   0,   0,   0,   1,   0,   0]
	issymmetric(A);
	L = extractL(A);
	[v,D] = eig(L);
	lambdas = diag(D)';
	eigenVectors = v;
	eigenValues = lambdas;
	FiedlerVector = v(:,2)'
	Lambda2 = lambdas(1,2) 
	printf("\n");
	for i=1:N
		fprintf(fid(i,1), "(%d, %f) \n", iter, v(i,2));
	endfor



	#centralized PI
	centralizedPIxs = (eye(N) - (1/N).*L)^50*L*centralizedPIxs;
	centralizedPIxs = centralizedPIxs./norm(centralizedPIxs);
	#centralizedPIFiedlerVector = centralizedPIxs'
	for i=1:N
		fprintf(fid(i,2), "(%d, %f) \n", iter, centralizedPIxs(i,1));
	endfor
	
	
	
	#distributed eigenvectors centrality
	Nk = sum(A,2);
	[m, index] = max(Nk);
	#update
	tempCs = zeros(N,1);
	for repeat=1:10
		for i=1:N
			tempCs(i,1) = A(i,:)*distrCs;
		endfor
		#tempCs = tempCs./norm(tempCs);
		tempCs = tempCs./(((A(index,:)*tempCs)/m)*N); #assumed unique over the graph and shared through VS
		distrCs = tempCs;
		tempCs = zeros(N,1);
	endfor
	#eigenvectorC = distrCs'



	#Bertrand distributed approach
	Nk = sum(A,2);
	[m, index] = max(Nk);
	#update
	distrXs = [0:9]'/10;
	distrXs = distrXs./norm(distrXs);
	tempXs = zeros(N,1);
	for i=1:N
		tempXs(i,1) = Nk(i,1)*distrXs(i,1)-A(i,:)*distrXs;
	endfor
	tempXs = tempXs./norm(tempXs);
	#tempXs = tempXs./(((norm(A(index,:).*tempXs'))/m)*N); #assumed unique over the graph and shared through VS
	distrXs = tempXs;
	tempXs = zeros(N,1);
	for repeat=1:100
		for i=1:N
			tempXs(i,1) = ((1-(Nk(i,1)/(3*N)))*distrXs(i,1))+(1/(3*N))*(A(i,:)*distrXs);
		endfor
		tempXs = tempXs./norm(tempXs);
		#tempXs = tempXs./(((norm(A(index,:).*tempXs'))/m)*N); #assumed unique over the graph and shared through VS
		distrXs = tempXs;
		tempXs = zeros(N,1);

	endfor
	bertrandFiedlerVector = (distrXs')#.*(FiedlerVector(1,1)/distrXs(1,1))
	bertrandDistributedLambda2s = (diag(L) - ((A*distrXs)./distrXs))'
	for i=1:N
		fprintf(fid(i,3), "(%d, %f) \n", iter, bertrandFiedlerVector(1,i));
	endfor
	
	
	
	#centralized DiLorenzo
	epsilon = 0.9*(2/N);
	for i=1:100
		centralizedDiLorenzo = ((eye(N) - epsilon.*L)-(1/N).*ones(N,N))*centralizedDiLorenzo;
		centralizedDiLorenzo = centralizedDiLorenzo./norm(centralizedDiLorenzo);
	endfor
	#centralizedDiLorenzo = centralizedDiLorenzo' 
	
	
	
	#DiLorenzo distributed approach
	Nk = sum(A,2);
	[m, index] = max(Nk);
	epsilon = 0.4*(2/N);
	oppA = abs(A-ones(N,N))-eye(N);
	#update
	distrBs = [0:9]'/10;
	distrBs = distrXs./norm(distrBs);
	tempBs = zeros(N,1);
	for repeat=1:100
		for i=1:N
			one = (((N-1)/N) - Nk(i,1)*epsilon)*distrBs(i,1);
			two = ((N*epsilon-1)/N).*A(i,:)*distrBs;
			three = (-1/N).*oppA(i,:)*distrBs;
			tempBs(i,1) = one + two + three;
			if i == 1
				one;
				two;
				three;
				oppA(i,:)*distrBs;
			endif
			#tempBs(i,1) = ((1-(Nk(i,1)/N))*distrBs(i,1))+(1/N)*(A(i,:)*distrBs);
		endfor
		tempBs = tempBs./norm(tempBs);
		#tempBs = tempBs./(((A(index,:)*tempBs)/m)*N); #assumed unique over the graph and shared through VS
		distrBs = tempBs;
		tempBs = zeros(N,1);
		#distrBs'
		#goOn = input("stop?");
	endfor
	diLorenzoFiedlerVector = (distrBs')#.*(FiedlerVector(1,1)/distrBs(1,1))
	diLorenzodDistributedLambda2s = (diag(L) - ((A*distrBs)./distrBs))'

	
		
	#Sahai distributed approach
	distrV = zeros(N,N);
	LNorm = zeros(N,N);
	for i=1:N
		LNorm(i,:) = L(i,:)./L(i,i);
	endfor
	if iter == 1
		for t=3:(Tmax+2)
			for n=1:N
				waveUs(t,n) = 2*waveUs(t-1,n) - waveUs(t-2,n) - c^2*(LNorm(n,:)*waveUs(t-1,:)');
			endfor
		endfor
	else
		waveUs(1:Tmax+1,:) = waveUs(2:Tmax+2,:);
		t = Tmax+2;
		for n=1:N
			waveUs(t,n) = 2*waveUs(t-1,n) - waveUs(t-2,n) - c^2*(LNorm(n,:)*waveUs(t-1,:)');
		endfor
	endif
	waveUs;
	#plot(1:(Tmax+2),waveUs(:,1),1:(Tmax+2),waveUs(:,2),1:(Tmax+2),waveUs(:,3),1:(Tmax+2),waveUs(:,4),1:(Tmax+2),waveUs(:,5),1:(Tmax+2),waveUs(:,6),1:(Tmax+2),waveUs(:,7),1:(Tmax+2),waveUs(:,8),1:(Tmax+2),waveUs(:,9),1:(Tmax+2),waveUs(:,10));
	#input("");
	Y = zeros(Tmax+2,N);
	for n=1:N
		Y(:,n) = fft(waveUs(:,n));
	endfor
	omegas = zeros(N,N);
	amp = zeros(N,N);
	phases = zeros(N,N);
	for i=1:N
		[peak_list, indexes] = findpeaks(abs(Y(1:(Tmax+2)/2,i)));
		omegas(i,1:length(indexes)) = 2*pi*(1/(Tmax+2)).*(indexes.-1);
		amp(i,1:length(peak_list)) = (1/(Tmax+2)).*peak_list;
		for j=1:length(indexes)
			phases(i,j) = arg(Y(indexes(j),i));
		endfor
	endfor
	waveDistributedOmega2s = omegas(:,2)'
	waveDistributedAmplitides = amp(:,2)'
	waveDistributedPhases = phases(:,2)'
	#waveDistributedLNormFiedler = 
	#waveDistributedLambda2s = (-1/c^2).*cos(omegas(:,2)').+(2/c^2)
	#waveDistributedLambda2s = (diag(L) - ((A*waveDistributedFiedlerVector')./waveDistributedFiedlerVector'))'
	[vNorm,DNorm] = eig(LNorm);
	LNormFiedlerVector = v(:,2)'
	LNormLambda2 = DNorm(2,2)

	#double check the trandform by rebuilding the signal
	whichnode = 1;
	rebuilt = zeros(1,Tmax+2);
	for i=1:Tmax+2
		a = abs(Y(i,whichnode));
		b = arg(Y(i,whichnode));
		for j=1:Tmax+2
			rebuilt(1,j) = rebuilt(1,j) + (1/(Tmax+2))*a*cos(2*pi*(1/(Tmax+2))*(i-1)*j+b);
		endfor
	endfor
	#subplot(2,1,1);
	#plot(1:Tmax+2,waveUs(:,whichnode));
	#subplot(2,1,2);
	#plot(1:Tmax+2,rebuilt(1,:));
	#input("");
	#plot(1:(Tmax+2)/2,abs(Y(1:(Tmax+2)/2,1)),1:(Tmax+2)/2,abs(Y(1:(Tmax+2)/2,2)),1:(Tmax+2)/2,abs(Y(1:(Tmax+2)/2,3)),1:(Tmax+2)/2,abs(Y(1:(Tmax+2)/2,4)),1:(Tmax+2)/2,abs(Y(1:(Tmax+2)/2,5)),1:(Tmax+2)/2,abs(Y(1:(Tmax+2)/2,6)),1:(Tmax+2)/2,abs(Y(1:(Tmax+2)/2,7)),1:(Tmax+2)/2,abs(Y(1:(Tmax+2)/2,8)),1:(Tmax+2)/2,abs(Y(1:(Tmax+2)/2,9)),1:(Tmax+2)/2,abs(Y(1:(Tmax+2)/2,10)));
	distrV;	



	goOn = input("proceed?");
until !isempty(goOn)
#until iter == 50

for i=1:N
	for j=1:3
		fputs(fid(i,j), "};");
		if j == 1 fputs(fid(i,j), "\\addlegendentry{actual} \n");
		elseif j == 2 fputs(fid(i,j), "\\addlegendentry{central} \n");
		else fputs(fid(i,j), "\\addlegendentry{distributed} \n");		
		endif
		fclose(fid(i,j));
	endfor
endfor





