function vPitch = estimatePitch(matFrames, fs)

	for i=1:rows(matFrames)
		flipped_frame =  fliplr(matFrames(i,:));
		acf(i,1:length(matFrames(i,:))*2-1) = conv(matFrames(i,:), flipped_frame);
	endfor

	zero_idx = length(matFrames(1,:))+1;
	start_idx = zero_idx + cast(1.0/400 * fs, "int16");
	end_idx = zero_idx + cast(1.0/80 * fs, "int16");
	for i=1:rows(acf)
		[max_val, max_index] = max(acf(i, start_idx:end_idx));
		true_max_idx = start_idx+max_index-1-zero_idx;
		f_period_len = double(true_max_idx)/double(fs);
		vPitch(i) = 1/f_period_len;
	endfor

endfunction
