#! /usr/bin/octave
#Ubuntu's interpreter /usr/bin/octave
#OSX's interpreter /usr/local/bin/octave

addpath("~/Dropbox/papers/2016/saso-carlo/");

	load ./matrices/A10s.mat; 
	#load ./matrices/A100s.mat; 
	#load ./matrices/A1000s.mat;

	#list parameters
	algorithmlist = [1];
	swarmsizelist = [1]; #[1,2,3];
	packetdroplist = [0.0]; #,0.01,0.05,0.1,0.2, 0.4,0.6,0.8];
	
	#scalar parameters
	maxiterations = 1000;
	errormodel = 1; #1: standard link failure with packet drop probability, 2: gentler error model, 3: ...
	
	#numoscillators = 3;
	numoscillators = 1;
	m = numoscillators;
	gam = 0.75;
	GAM = ones(1,m);
	
	function retval = EF(x) #dynamics of the uncoupled system (gam -> 0), (assumed chaotic ?)
	  retval = 0;
	  if (isvector (x))
	    #retval = ones(rows(x),columns(x))./x;
	    if length(x) == 1
	    	a2 = 1.3;
	    	retval = (-1.0*a2*sign(x)*(1 - x))-x;
	    endif
	    if length(x) == 3 
	    	s = 10;
	    	r = 28;
	    	b = 8/3;
	    	step = 100;
	    	retval = zeros(1,3);
	    	retval(1) = (s*(x(2)-x(1)))/step;
	    	retval(2) = (r*x(1) - x(1)*x(3) - x(2))/step;
	    	retval(3) = (x(1)*x(2) - b*x(3))/step;
	    endif 
	  else
	    error ("EF: expecting vector argument");
	  endif
	endfunction
	
	function retval = AITCH(x)
	  retval = 0;
	  if (isvector (x))
	    retval = mean(x);
	  else
	    error ("AITCH: expecting vector argument");
	  endif
	endfunction
	

	fileids = zeros(m,1);
	for algorithm=algorithmlist
		switch (algorithm)
	  	case 1 name = "sorrentino";
		#
		#
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
				#
				#
				endswitch
				for o=1:numoscillators
					fileids(o,1) = fopen(strcat("./data/oscillator-",int2str(o),"-errorvssteps-",name,"-",int2str(swarmsize),"-",num2str(packetdrop),"-",errormodelname,".dat"), "w");
				endfor
				avgtime = 0.0;
				printf(cstrcat("\n", name, " swarm size=", int2str(swarmsize), " drop p=", num2str(packetdrop), " (error model ", num2str(errormodel), ")\n"));
				for sample=[1] #[1,..,length(currentCell)]
					tic ();
					currentA = currentCell{sample};
					
					#make the adjacency matrix fully connected
					#currentA = ones(rows(currentA), rows(currentA));
					#currentA = currentA - eye(rows(currentA))
					
					#from Sorrentino, 2015
					currentA = [0, 0, 0, 1, 0; 
						    0, 0, 0, 1, 0; 
						    0, 0, 0, 1, 1; 
						    1, 1, 1, 0, 0;
						    0, 0, 1, 0, 0];
						    
					#compare
					#currentA = [1, 1, 1, 1; 1, 1, 1, 1;  1, 1, 1, 1; 1, 1, 1, 1]
					
					
					currentK = rows(currentA);
					
					for o=1:numoscillators
						fprintf(fileids(o,1), "i ");
						for node=1:9 #rows(currentA)-1
							fprintf(fileids(o,1), "node%d ", node);
						endfor
						fprintf(fileids(o,1), "node%d\n", 10);
					endfor
					#
					#
					
					#identical starting point
					Xs = 0.15*ones(rows(currentA), m);
					#random starting points
					Xs = 0.15*rand(rows(currentA), m);
					
					###
					for iteration=1:maxiterations
						###
						#
						#
						printf("%d\r", iteration);
						
	#input("");
						
						
						if iteration == 1200
							currentA = [0, 0, 0, 1; 
						    		    0, 0, 0, 1; 
						    		    0, 0, 0, 0; 
						    		    1, 1, 0, 0]
						endif
						
						if iteration == 1400
							currentA = [0, 0, 0, 1; 
						    		    0, 0, 0, 1; 
						    		    0, 0, 0, 1; 
						    		    1, 1, 1, 0]
						endif
						
						###
						#run algorithm
					  	switch (algorithm)
						case 1
							###
							#
							# EVOLVE Xs
							
							rs = zeros(rows(currentA),1);
							droppedpacket = zeros(rows(currentA),1);
							for node=1:rows(currentA)
								for neigh=1:columns(currentA)
								if currentA(node,neigh) != 0.0
								packetdrop = rand();
									if packetdrop >= 0.0 #&& node != neigh #ignore the self-loop?
										rs(node,1) += currentA(node,neigh).*AITCH(Xs(neigh,:));
									else
										droppedpacket(node,1) += 1;
									endif
								endif
								endfor
							endfor
							rs
							
							#sigmas = ones(rows(currentA),1) ./ sum(currentA, 2);
							# recompute sigma to account for packet drops, i.e. DEAD LINKS
							sigmas = ones(rows(currentA),1) ./ (sum(currentA, 2)-droppedpacket);
							
							sigmas(isinf(sigmas)) = 0; 
							sigmas
							
							for node=1:rows(currentA)
							
								#uncoupled dynamics
								#XsPoint = EF(Xs(node,:));
								
								if sigmas(node,1) != 0
								XsPoint = EF(Xs(node,:)) + gam.*GAM .*( sigmas(node,1) * rs(node,1) - AITCH(Xs(node,:)) );
								else
								XsPoint = EF(Xs(node,:));
								endif
								Xs(node,:) = Xs(node,:) + XsPoint;
							endfor
							Xs
							
								
							#
						endswitch
						###
						
						#write files
						if mod(iteration, 10) == 1
							for o=1:numoscillators
								fprintf(fileids(o,1), "%f ", iteration);
								for node=1:9 #rows(currentA)-1
									if node <= rows(currentA)
									fprintf(fileids(o,1), "%f ", Xs(node,o));
									else
									fprintf(fileids(o,1), "%f ", 0.0);
									endif
								endfor
								if 10 == rows(currentA)
								fprintf(fileids(o,1), "%f\n", Xs(rows(currentA),o));
								else
									fprintf(fileids(o,1), "%f\n", 0.0);
								endif
							endfor
						endif
						
					endfor
					avgtime = (avgtime*(sample-1) + toc())/sample;
				endfor
				###
				#close files
				for o=1:numoscillators
					fclose(fileids(o,1));
				endfor
				###
			endfor
		endfor
	endfor
