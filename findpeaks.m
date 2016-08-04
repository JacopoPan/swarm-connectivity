function [peak_list, indexes] = findpeaks(array)
	if length(array) == 0
		printf("using findpeaks on empty array\n");
	elseif length(array) == 1
		indexes = 1;
		peak_list = array(1);
	else
		counter = 0;
		if array(1) > array(2)
			counter++;
		endif
		for i=2:length(array)-1
			if array(i) > array(i-1) && array(i) > array(i+1)
				counter++;
			endif
		endfor
		if array(length(array)) > array(length(array)-1)
			counter++;
		endif
	
		indexes = zeros(1,counter);
		peak_list = zeros(1,counter);
		index = 1;
		if array(1) > array(2)
			indexes(index) = 1;
			peak_list(index) = array(1);
			index = index+1;
		endif
		for i=2:length(array)-1
			if array(i) > array(i-1) && array(i) > array(i+1)
				indexes(index) = i;
				peak_list(index) = array(i);
				index = index+1;
			endif
		endfor
		if array(length(array)) > array(length(array)-1)
			indexes(index) = length(array);
			peak_list(index) = array(length(array));
		endif
	endif



	#optional removal of non-relevant peaks



endfunction
