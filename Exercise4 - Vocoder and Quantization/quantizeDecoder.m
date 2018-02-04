function qx = quantizeDecoder(idx,nBits,xOffset,xCenter)
	len = (xCenter+xOffset) - (xCenter-xOffset);
	stepSize = len / (2^nBits - 1);
	q_levels = (xCenter-xOffset) : stepSize : (xCenter+xOffset);
	for i=1:length(idx)
		qx(i) = q_levels(idx(i));
	endfor
endfunction
