function newX = consensusIteration (A, currentX)
	K = rows(A);
	gam = 0.8*(1/(K-1));
	D = diag(sum(A,2));
	matrix = eye(K) - gam.*D + gam.*A;
	newX = matrix*currentX;
endfunction