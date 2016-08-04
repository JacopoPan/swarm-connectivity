function omega = preComputeOmega(A,maxiter)
	K = rows(A);
	Us = zeros(2*maxiter,K); 
	Us(2,:) = unifrnd(0.0,1.0,[1,K]); 
	Us(1,:) = Us(2,:);
	
	L = extractL(A);
	LNorm = zeros(K,K);
	for k=1:K
		LNorm(k,:) = L(k,:)./L(k,k);
	endfor
	
	for t = 3:rows(Us)
		for k=1:K
			Us(t,k) = 2*Us(t-1,k) - Us(t-2,k) - 0.8^2*(LNorm(k,:)*Us(t-1,:)');
		endfor
	endfor
	
	Y = zeros(2*maxiter,K);
	for k=1:K
		Y(:,k) = fft(Us(:,k));
	endfor
	
	omegas = zeros(K,K);
	for k=1:K
		[peak_list, indexes] = findpeaks(abs(Y(1:(2*maxiter)/2,k)));
		omegas(k,1:length(indexes)) = 2*pi*(1/(2*maxiter)).*(indexes.-1);
	endfor
	omega = mode(omegas(:,2)');
endfunction
