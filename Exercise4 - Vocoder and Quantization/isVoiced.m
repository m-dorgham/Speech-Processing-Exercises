function seg_voiced = isVoiced(seg)
	exactZeros = find(seg==0);
	newSeg = seg(1:end-1).*seg(2:end);
	signChange = find(newSeg<0);
	zcCount = length(exactZeros) + length(signChange);
	zcRate = zcCount / length(seg);
	segPower = sqrt(computePower(seg));
	if (zcRate<.3 && segPower>.01)
		seg_voiced = 1;
	else
		seg_voiced = 0;
	endif
endfunction
