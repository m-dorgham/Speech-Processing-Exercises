[x,fs] = audioread ("female8khz.wav");
t= linspace(0, length(x)/fs, length(x));

frame_length = 0.032;
frame_shift = 0.008;

[m_frames, v_time_frame] = my_windowing(x, fs, frame_length, frame_shift);

for i=1:rows(m_frames)
	x_energy(i) = sqrt(computePower(m_frames(i,:)));
endfor

figure
plot(t, x)
hold on
plot(v_time_frame/fs, x_energy, 'LineWidth',2)
hold off
axis tight


for i=1:rows(m_frames)
	seg_voiced(i) = isVoiced(m_frames(i,:));
endfor

figure
plot(t, x)
hold on
plot(v_time_frame/fs, seg_voiced.*0.1, 'LineWidth',2)
hold off
axis tight


vPitch = estimatePitch(m_frames, fs);
[m_stft, v_freq, v_time] = my_stft(x, fs, frame_length, frame_shift, hann(frame_length*fs, 'periodic'));

figure
imagesc(v_time/fs, v_freq, 10*log10(max(abs(m_stft).^2, 10^(-15))))
axis xy
hold on
plot(v_time/fs, vPitch, 'r', 'LineWidth',2)
hold off



matLPC = getLPC(m_frames, 15);

%[H, freq_ax] = freqz(1,[1;matLPC(:,30)],columns(m_frames),'whole',fs);
%S = fft(m_frames(30,:)'.*hann(columns(m_frames), 'periodic'));

%figure
%plot(freq_ax, 20*log10(abs(H)))
%hold on
%plot(freq_ax, 20*log10(abs(S)))
%hold off
%xlabel('Frequency (Hz)')
%ylabel('Magnitude (dB)')
%axis tight

