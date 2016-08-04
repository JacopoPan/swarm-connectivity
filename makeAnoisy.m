function nA = makeAnoisy (A,p)
	N = rows(A);
	ra = rand(N,N);
	#nA = A - A.*(rand(N,N)<p);
	nA = A.*(rand(N,N)>p);
		
	#nA = A;
	#[is,js] = find(nA);
	#for n=1:length(is)
	#	if unifrnd(0,1) < p
	#		nA(is(n),js(n)) = 0;
	#	endif
	#endfor
	
	#for i=1:rows(nA)
	#	for j=i:columns(nA)
	#		if nA(i,j) == 1 && unifrnd(0,1) < p
	#			nA(i,j) = 0;
	#		endif
	#	endfor
	#endfor
endfunction