[M_1,fs_1] = audioread ("/home/dorgham/Documents/AudioFiles/speech1.wav");
[M_2,fs_2] = audioread ("/home/dorgham/Documents/AudioFiles/speech2.wav");


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


[m_frames_1, v_time_frame_1] = my_windowing(M_1, fs_1, 0.032, 0.016);
[m_frames_2, v_time_frame_2] = my_windowing(M_2, fs_2, 0.032, 0.016);


for i=1:length(m_frames_1(:,1))
	flipped_frame =  fliplr(m_frames_1(i,:));
	acf_1(i,1:length(m_frames_1(i,:))*2-1) = conv(m_frames_1(i,:), flipped_frame);
endfor


for i=1:length(m_frames_2(:,1))
	flipped_frame =  fliplr(m_frames_2(i,:));
	acf_2(i,1:length(m_frames_2(i,:))*2-1) = conv(m_frames_2(i,:), flipped_frame);
endfor


zero_idx = length(m_frames_1(1,:))+1;
start_idx = zero_idx + cast(1.0/400 * fs_1, "int16");
end_idx = zero_idx + cast(1.0/80 * fs_1, "int16");
for i=1:length(acf_1(:,1))
	[max_val, max_index] = max(acf_1(i, start_idx:end_idx));
	true_max_idx = start_idx+max_index-1-zero_idx;
	f_period_len = double(true_max_idx)/double(fs_1);
	f_freqs_1(i) = 1/f_period_len;
	X = ["Fundemental frequency is ", num2str(f_freqs_1(i)), " Hz. At lag of ", num2str(f_periods_1(i)*1000), " ms."];
endfor

zero_idx = length(m_frames_2(1,:))+1;
start_idx = zero_idx + cast(1.0/400 * fs_2, "int16");
end_idx = zero_idx + cast(1.0/80 * fs_2, "int16");
for i=1:length(acf_2(:,1))
	[max_val, max_index] = max(acf_2(i, start_idx:end_idx));
	true_max_idx = start_idx+max_index-1-zero_idx;
	f_period_len = double(true_max_idx)/double(fs_2);
	f_freqs_2(i) = 1/f_period_len;
	X = ["Fundemental frequency is ", num2str(f_freqs_2(i)), " Hz. At lag of ", num2str(f_periods_2(i)*1000), " ms."];
endfor


% plot the estimated pitches

t1= linspace(0, length(M_1)/fs_1, length(M_1));
t2= linspace(0, length(M_2)/fs_2, length(M_2));
x1= v_time_frame_1/fs_1;
x2= v_time_frame_2/fs_2;

figure
subplot(2,1,1)       % add first plot in 2 x 1 grid
plot(t1,M_1)
title('Speech Signal 1')

subplot(2,1,2)       % add second plot in 2 x 1 grid
plot(x1,f_freqs_1)       
title('Estimated fundamental frequency')


figure
subplot(2,1,1)       % add first plot in 2 x 1 grid
plot(t2,M_2)
title('Speech Signal 2')

subplot(2,1,2)       % add second plot in 2 x 1 grid
plot(x2,f_freqs_2)       
title('Estimated fundamental frequency')

