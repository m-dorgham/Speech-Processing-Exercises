function [m_stft, v_freq, v_time] = my_stft(v_signal, sampling_rate, frame_length, frame_shift, v_analysis_window)
	m_frames_width = frame_length*sampling_rate;
	m_frames_height= length(v_signal)/(frame_shift*sampling_rate);
	v_freq = sampling_rate*(0:m_frames_width/2)/m_frames_width;
	i = 1;
	for j = 1:frame_shift*sampling_rate:length(v_signal)
		v_time(i)=j+m_frames_width-(m_frames_width/2)-1;
		if(j+m_frames_width-1 < length(v_signal))
			m_frames(i, 1:m_frames_width) = v_signal(j:j+m_frames_width-1);
		else
			remaining_elems = length(v_signal)-j;
			m_frames(i, 1:remaining_elems) = v_signal(j:j+remaining_elems-1);
			m_frames(i, remaining_elems+1:m_frames_width) = zeros(1, m_frames_width-remaining_elems);
		endif
		%compute windowed signal
		x = m_frames(i,:)'.*v_analysis_window;
		%compute fft
		Y = fft(x);
		m_stft(:,i) = Y(1:length(Y)/2+1);
		i++;
	endfor
endfunction
