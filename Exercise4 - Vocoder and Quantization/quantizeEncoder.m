function idx = quantizeEncoder(x,nBits,xOffset,xCenter)
	len = (xCenter+xOffset) - (xCenter-xOffset);
	stepSize = len / (2^nBits - 1);
	q_levels = (xCenter-xOffset) : stepSize : (xCenter+xOffset);
	for i=1:length(x)
		[m,j] = min(abs(x(i)-q_levels));
		idx(i) = j;
	endfor
endfunction
