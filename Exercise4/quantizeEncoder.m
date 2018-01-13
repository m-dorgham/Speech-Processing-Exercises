function idx = quantizeEncoder(x,nBits,xMax,xCenter)
	len = xMax - (2*xCenter-xMax);
	stepSize = len / (2^nBits - 1);
	q_levels = (2*xCenter-xMax) : stepSize : xMax;
	for i=1:length(x)
		[m,j] = min(abs(x(i).-q_levels));
		idx(i) = j;
	endfor
endfunction
