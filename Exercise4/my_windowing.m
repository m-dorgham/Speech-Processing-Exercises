function [m_frames, v_time_frame] = my_windowing(v_signal, sampling_rate, frame_length, frame_shift)
	m_frames_width = frame_length*sampling_rate;
	m_frames_height= length(v_signal)/(frame_shift*sampling_rate);
	i = 1;
	for j = 1:frame_shift*sampling_rate:length(v_signal)
		v_time_frame(i)=j+m_frames_width-(m_frames_width/2)-1;
		if(j+m_frames_width-1 < length(v_signal))
			m_frames(i, 1:m_frames_width) = v_signal(j:j+m_frames_width-1);
		else
			remaining_elems = length(v_signal)-j;
			m_frames(i, 1:remaining_elems) = v_signal(j:j+remaining_elems-1);
			m_frames(i, remaining_elems+1:m_frames_width) = zeros(1, m_frames_width-remaining_elems);
		endif
		i++;
	endfor

endfunction


