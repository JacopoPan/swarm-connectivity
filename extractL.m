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