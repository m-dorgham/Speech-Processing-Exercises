function s_power = computePower(x)
	s_power = sum(abs(x).^2)/length(x);
endfunction
