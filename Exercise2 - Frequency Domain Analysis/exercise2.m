close all
pkg load signal

[M_1,fs_1] = audioread ("Audio/speech1.wav");
[M_2,fs_2] = audioread ("Audio/phone.wav");


%subsection 1
frame_length=0.032;
frame_shift=0.008;
[m_stft_1, v_freq_1, v_time_1] = my_stft(M_1, fs_1, frame_length, frame_shift, hann(frame_length*fs_1, 'periodic'));
[m_stft_2, v_freq_2, v_time_2] = my_stft(M_2, fs_2, frame_length, frame_shift, hann(frame_length*fs_2, 'periodic'));

%subsection 2
[m_frames_1, v_time_frame_1] = my_windowing(M_1, fs_1, 0.032, 0.016);
for i=1:rows(m_frames_1)
	flipped_frame =  fliplr(m_frames_1(i,:));
	acf_1(i,1:length(m_frames_1(i,:))*2-1) = conv(m_frames_1(i,:), flipped_frame);
endfor
zero_idx = columns(m_frames_1)+1;
start_idx = zero_idx + cast(1.0/400 * fs_1, "int16");
end_idx = zero_idx + cast(1.0/80 * fs_1, "int16");
for i=1:rows(acf_1)
	[max_val, max_index] = max(acf_1(i, start_idx:end_idx));
	true_max_idx = start_idx+max_index-1-zero_idx;
	f_period_len = double(true_max_idx)/double(fs_1);
	f_freqs_1(i) = 1/f_period_len;
endfor
x1= v_time_frame_1/fs_1;


figure
imagesc(v_time_1/fs_1, v_freq_1, 10*log10(max(abs(m_stft_1).^2, 10^(-15))))
axis xy
xlabel('time (s)')
ylabel('Frequency (Hz)')
hold on
plot(x1,f_freqs_1, 'r', 'LineWidth',2)
hold off
title('Spectrogram of the speech signal with estimated f0')
drawnow

figure
imagesc(v_time_2/fs_2, v_freq_2, 10*log10(max(abs(m_stft_2).^2, 10^(-15))))
axis xy
xlabel('time (s)')
ylabel('Frequency (Hz)')
title('Spectrogram of the phone signal')
drawnow

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

