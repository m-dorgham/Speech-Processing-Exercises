function matLPC = getLPC(matFrames, M)
	
	for i=1:rows(matFrames)
		win_frame = matFrames(i,:)'.*hann(columns(matFrames), 'periodic');
		%compute the autocorrelation vector
		phi_s = xcorr(win_frame);
		phi_s = phi_s(round(length(phi_s)/2):round(length(phi_s)/2)+M);
		%compute the toeplitz matrix
		R_s = toeplitz(phi_s(1:end-1));
		phi_s = phi_s(2:end);
		%compute LP coefficients
		a = -(R_s\phi_s);
		matLPC(:,i) = a';
	endfor

endfunction
