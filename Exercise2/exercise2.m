[M_1,fs_1] = audioread ("/home/dorgham/Documents/AudioFiles2/speech1.wav");
[M_2,fs_2] = audioread ("/home/dorgham/Documents/AudioFiles2/phone.wav");


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

function v_signal = my_inverse_stft_2(m_stft, fs, frame_length, frame_shift, v_synthesis_window)
	frame_len_samples = frame_length*fs;
	frame_shft_samples = frame_shift*fs;
	cols_n = length(m_stft(1,:));
	v_sig_len = frame_len_samples+(cols_n-1)*frame_shft_samples;
	v_signal = zeros(1, v_sig_len);

	index=1;
	for i=1:cols_n
		%complete the stft with its mirror conjugate
		v_stft = m_stft(:, i);
		v_complete_stft = [v_stft; conj(v_stft(end-1:-1:2))];
		%compute inverse fft for the current column
		v_istft = ifft(v_complete_stft);
		%overlap-add method
		v_signal(index:index+frame_len_samples-1) = v_signal(index:index+frame_len_samples-1) + real(v_istft)'.*v_synthesis_window';
		index = index + frame_shft_samples;
	endfor
endfunction

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

%subsection 1
frame_length=0.032;
frame_shift=0.008;
[m_stft_1, v_freq_1, v_time_1] = my_stft(M_1, fs_1, frame_length, frame_shift, hann(frame_length*fs_1, 'periodic'));
[m_stft_2, v_freq_2, v_time_2] = my_stft(M_2, fs_2, frame_length, frame_shift, hann(frame_length*fs_2, 'periodic'));

%subsection 2
[m_frames_1, v_time_frame_1] = my_windowing(M_1, fs_1, 0.032, 0.016);
for i=1:length(m_frames_1(:,1))
	flipped_frame =  fliplr(m_frames_1(i,:));
	acf_1(i,1:length(m_frames_1(i,:))*2-1) = conv(m_frames_1(i,:), flipped_frame);
endfor
zero_idx = length(m_frames_1(1,:))+1;
start_idx = zero_idx + cast(1.0/400 * fs_1, "int16");
end_idx = zero_idx + cast(1.0/80 * fs_1, "int16");
for i=1:length(acf_1(:,1))
	[max_val, max_index] = max(acf_1(i, start_idx:end_idx));
	true_max_idx = start_idx+max_index-1-zero_idx;
	f_period_len = double(true_max_idx)/double(fs_1);
	f_freqs_1(i) = 1/f_period_len;
endfor
x1= v_time_frame_1/fs_1;


figure
imagesc(v_time_1/fs_1, v_freq_1, 10*log10(max(abs(m_stft_1).^2, 10^(-15))))
axis xy
hold on
plot(x1,f_freqs_1, 'r', 'LineWidth',2)
hold off

figure
imagesc(v_time_2/fs_2, v_freq_2, 10*log10(max(abs(m_stft_2).^2, 10^(-15))))
axis xy

%subsection 3
v_test_signal = ones(2048, 1);
fs_test = 16000;

[m_stft_test, v_freq_test, v_time_test] = my_stft(v_test_signal, fs_test, 0.032, 0.016, sqrt(hann(0.032*fs_test, 'periodic')));

v_synth_signal = my_inverse_stft(m_stft_test, fs_test, 0.032, 0.016, sqrt(hann(0.032*fs_test, 'periodic')));

%plot the synthesized signal
t1= linspace(0, length(v_test_signal)/fs_test, length(v_test_signal));
t2 = linspace(0, length(v_synth_signal)/fs_test, length(v_synth_signal));

figure

ax1 = subplot(2,1,1);
plot(ax1, t1, v_test_signal)
title('Original Signal')
axis(ax1, [0 length(v_synth_signal)/fs_test -inf inf])

ax2 = subplot(2,1,2);
plot(ax2, t2, v_synth_signal)
title('Synthesized Signal')
axis(ax2,[0 length(v_synth_signal)/fs_test -inf inf])

