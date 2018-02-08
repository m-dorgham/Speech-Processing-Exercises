close all
pkg load signal

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
title('The energy of the speech signal')
drawnow

for i=1:rows(m_frames)
	seg_voiced(i) = isVoiced(m_frames(i,:));
endfor

figure
plot(t, x)
hold on
plot(v_time_frame/fs, seg_voiced.*0.1, 'LineWidth',2)
hold off
axis tight
title('Voiced and Unvoiced areas in the signal')
drawnow

vPitch = estimatePitch(m_frames, fs);
[m_stft, v_freq, v_time] = my_stft(x, fs, frame_length, frame_shift, hann(frame_length*fs, 'periodic'));

figure
imagesc(v_time/fs, v_freq, 10*log10(max(abs(m_stft).^2, 10^(-15))))
axis xy
xlabel('time (s)')
ylabel('Frequency (Hz)')
colorbar
hold on
plot(v_time/fs, vPitch, 'r', 'LineWidth',2)
hold off
title('Spectrogram of the signal with the estimated fundamental frequency')
drawnow


matLPC = getLPC(m_frames, 12);

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

%%%% 3.1
nSamples = 18000;
nSamplesPerSecond = nSamples / 2.25;
jump = nSamplesPerSecond / 100;
e_unvoiced = randn(nSamples,1);
e_voiced = zeros(nSamples, 1);
e_voiced(1:jump:end) = 1;

[m_frames_v, v_time_frame_v] = my_windowing(e_voiced, nSamplesPerSecond, frame_shift, frame_shift);
[m_frames_u, v_time_frame_u] = my_windowing(e_unvoiced, nSamplesPerSecond, frame_shift, frame_shift);

filterState = [];
signal_out_v = [];
signal_out_u = [];
for i=1:rows(m_frames_v)
	[segmentOut_v,filterState] = filterAdaptively(1,[1;matLPC(:,i)],m_frames_v(i,:),filterState);
	signal_out_v = [signal_out_v; segmentOut_v];
endfor

filterState = [];
for i=1:rows(m_frames_u)
	[segmentOut_u,filterState] = filterAdaptively(1,[1;matLPC(:,i)],m_frames_u(i,:),filterState);
	signal_out_u = [signal_out_u; segmentOut_u];
endfor

%play the sounds
playObj = audioplayer(signal_out_v,nSamplesPerSecond);
%playblocking(playObj);
playObj = audioplayer(signal_out_u,nSamplesPerSecond);
%playblocking(playObj);


[m_stft_v, v_freq_v, v_time_v] = my_stft(signal_out_v, nSamplesPerSecond, frame_length, frame_shift, hann(frame_length*nSamplesPerSecond, 'periodic'));
[m_stft_u, v_freq_u, v_time_u] = my_stft(signal_out_u, nSamplesPerSecond, frame_length, frame_shift, hann(frame_length*nSamplesPerSecond, 'periodic'));
figure
set(gcf, 'Position', [100, 200, 1100, 400])
subplot(1,2,1);
imagesc(v_time_v/nSamplesPerSecond, v_freq_v, 10*log10(max(abs(m_stft_v).^2, 10^(-15))))
axis xy
xlabel('time (s)')
ylabel('Frequency (Hz)')
colorbar
title('Spectrogram of generated voiced signal')
subplot(1,2,2);
imagesc(v_time_u/nSamplesPerSecond, v_freq_u, 10*log10(max(abs(m_stft_u).^2, 10^(-15))))
axis xy
xlabel('time (s)')
ylabel('Frequency (Hz)')
colorbar
title('Spectrogram of generated unvoiced signal')
drawnow

%%3.2
filterState = [];
signal_out = [];
for i=1:rows(m_frames_v)
	if seg_voiced(i)==1
		sel_frame = m_frames_v(i,:);
	else 
		sel_frame = m_frames_u(i,:);
	end
	[segmentOut, filterState] = filterAdaptively(1, [1;matLPC(:,i)], sel_frame, filterState);
	signal_out = [signal_out; segmentOut];
endfor

playObj = audioplayer(signal_out, nSamplesPerSecond);
%playblocking(playObj);
%plot the spectrogram
[m_stft_c, v_freq_c, v_time_c] = my_stft(signal_out, nSamplesPerSecond, frame_length, frame_shift, hann(frame_length*nSamplesPerSecond, 'periodic'));
figure
imagesc(v_time_c/nSamplesPerSecond, v_freq_c, 10*log10(max(abs(m_stft_c).^2, 10^(-15))))
axis xy
xlabel('time (s)')
ylabel('Frequency (Hz)')
title('Incorporating voiced/unvoiced information')
drawnow

%%3.3 Adjusting the amplitude
for i=1:rows(m_frames_v)
	frame_energy_v = max(sqrt(computePower(m_frames_v(i,:))),10^(-15));
	frame_energy_u = max(sqrt(computePower(m_frames_u(i,:))),10^(-15));
	gain_v = x_energy(i)/frame_energy_v;
	gain_u = x_energy(i)/frame_energy_u;
	m_frames_v(i,:) = m_frames_v(i,:).*gain_v;
	m_frames_u(i,:) = m_frames_u(i,:).*gain_u;
endfor

filterState = [];
signal_out = [];
for i=1:rows(m_frames_v)
	if seg_voiced(i)==1
		sel_frame = m_frames_v(i,:);
	else 
		sel_frame = m_frames_u(i,:);
	end
	[segmentOut, filterState] = filterAdaptively(1, [1;matLPC(:,i)], sel_frame, filterState);
	signal_out = [signal_out; segmentOut];
endfor

playObj = audioplayer(signal_out, nSamplesPerSecond);
%playblocking(playObj);

%plot the spectrogram
[m_stft_g, v_freq_g, v_time_g] = my_stft(signal_out, nSamplesPerSecond, frame_length, frame_shift, hann(frame_length*nSamplesPerSecond, 'periodic'));
figure
imagesc(v_time_g/nSamplesPerSecond, v_freq_g, 10*log10(max(abs(m_stft_g).^2, 10^(-15))))
axis xy
xlabel('time (s)')
ylabel('Frequency (Hz)')
title('Adjusting the generated signal power')
drawnow

%3.4 involve fundamental frequency
f_freq_samples = nSamplesPerSecond ./ vPitch;
frame_len_samples = columns(m_frames_v);

e_voiced = zeros(nSamples, 1);
counter = 0;
for i=1:length(e_voiced)
	counter = counter+1;
	if counter>=f_freq_samples(ceil(i/frame_len_samples))
		counter=0;
		e_voiced(i)=1;
	end
end

[m_frames_f, v_time_frame_f] = my_windowing(e_voiced, nSamplesPerSecond, frame_shift, frame_shift);

% Adjust the amplitude
for i=1:rows(m_frames_f)
	frame_energy_f = max(sqrt(computePower(m_frames_f(i,:))),10^(-15));
	gain_f = x_energy(i)/frame_energy_f;
	m_frames_f(i,:) = m_frames_f(i,:).*gain_f;
endfor

%apply adaptive filtering
filterState = [];
signal_out_f = [];
for i=1:rows(m_frames_f)
	if seg_voiced(i)==1
		sel_frame = m_frames_f(i,:);
	else 
		sel_frame = m_frames_u(i,:);
	end
	[segmentOut, filterState] = filterAdaptively(1, [1;matLPC(:,i)], sel_frame, filterState);
	signal_out_f = [signal_out_f; segmentOut];
endfor

playObj = audioplayer(signal_out_f, nSamplesPerSecond);
playblocking(playObj);

[m_frames_o, v_time_frame_o] = my_windowing(signal_out_f, nSamplesPerSecond, frame_length, frame_shift);
vPitch_f = estimatePitch(m_frames_o, nSamplesPerSecond);

%plot the spectrogram
[m_stft_f, v_freq_f, v_time_f] = my_stft(signal_out_f, nSamplesPerSecond, frame_length, frame_shift, hann(frame_length*nSamplesPerSecond, 'periodic'));
figure
imagesc(v_time_f/nSamplesPerSecond, v_freq_f, 10*log10(max(abs(m_stft_f).^2, 10^(-15))))
axis xy
xlabel('time (s)')
ylabel('Frequency (Hz)')
hold on
plot(v_time_f/nSamplesPerSecond, vPitch_f, 'r', 'LineWidth',2)
hold off
title('Incorporating fundamental frequency of original signal')
drawnow

%Quantization test with ramp function
rs=-5:0.01:5;
xMax=3;
xCenter=0;
nBits=2;

enc_rs = quantizeEncoder(rs, nBits, xMax, xCenter);
q_rs = quantizeDecoder(enc_rs, nBits, xMax, xCenter);

figure
plot(rs,'DisplayName','ramp function')
hold on
plot(q_rs,'DisplayName','quantized ramp')
hold off
axis tight
title('Reconstructing quantized ramp function')
legend('show')
drawnow

%Quantizing fundamental frequency
xCenter = round(mean(vPitch));
xMax = round(max(vPitch) - mean(vPitch));
nBits = 6;

enc_vPitch = quantizeEncoder(vPitch, nBits, xMax, xCenter);
q_vPitch = quantizeDecoder(enc_vPitch, nBits, xMax, xCenter);

%resynthesize the signal
f_freq_samples = nSamplesPerSecond ./ q_vPitch;
frame_len_samples = columns(m_frames_v);

e_voiced = zeros(nSamples, 1);
counter = 0;
for i=1:length(e_voiced)
	counter = counter+1;
	if counter>=f_freq_samples(ceil(i/frame_len_samples))
		counter=0;
		e_voiced(i)=1;
	end
end

[m_frames_f, v_time_frame_f] = my_windowing(e_voiced, nSamplesPerSecond, frame_shift, frame_shift);

% Adjust the amplitude
for i=1:rows(m_frames_f)
	frame_energy_f = max(sqrt(computePower(m_frames_f(i,:))),10^(-15));
	gain_f = x_energy(i)/frame_energy_f;
	m_frames_f(i,:) = m_frames_f(i,:).*gain_f;
endfor

%apply adaptive filtering
filterState = [];
signal_out_f = [];
for i=1:rows(m_frames_f)
	if seg_voiced(i)==1
		sel_frame = m_frames_f(i,:);
	else 
		sel_frame = m_frames_u(i,:);
	end
	[segmentOut, filterState] = filterAdaptively(1, [1;matLPC(:,i)], sel_frame, filterState);
	signal_out_f = [signal_out_f; segmentOut];
endfor

playObj = audioplayer(signal_out_f, nSamplesPerSecond);
playblocking(playObj);

pitch_combined=[vPitch; q_vPitch]';
figure
hist(pitch_combined,500)
title('Histogram of quantized and unquantized fundamental frequency')

%Quantizing the signal energy
figure
hist(x_energy)
title('Histogram of signal energy per frames')

xCenter = (max(x_energy)-min(x_energy))/2;
xMax = max(x_energy) - xCenter;
nBits = 3;

enc_energy = quantizeEncoder(x_energy, nBits, xMax, xCenter);
q_energy = quantizeDecoder(enc_energy, nBits, xMax, xCenter);

%resynthesize the signal
f_freq_samples = nSamplesPerSecond ./ q_vPitch;
frame_len_samples = columns(m_frames_v);

e_voiced = zeros(nSamples, 1);
counter = 0;
for i=1:length(e_voiced)
	counter = counter+1;
	if counter>=f_freq_samples(ceil(i/frame_len_samples))
		counter=0;
		e_voiced(i)=1;
	end
end

[m_frames_f, v_time_frame_f] = my_windowing(e_voiced, nSamplesPerSecond, frame_shift, frame_shift);

% Adjust the amplitude
for i=1:rows(m_frames_f)
	frame_energy_f = max(sqrt(computePower(m_frames_f(i,:))),10^(-15));
	gain_f = q_energy(i)/frame_energy_f;
	m_frames_f(i,:) = m_frames_f(i,:).*gain_f;
endfor

%apply adaptive filtering
filterState = [];
signal_out_f = [];
for i=1:rows(m_frames_f)
	if seg_voiced(i)==1
		sel_frame = m_frames_f(i,:);
	else 
		sel_frame = m_frames_u(i,:);
	end
	[segmentOut, filterState] = filterAdaptively(1, [1;matLPC(:,i)], sel_frame, filterState);
	signal_out_f = [signal_out_f; segmentOut];
endfor

playObj = audioplayer(signal_out_f, nSamplesPerSecond);
playblocking(playObj);


log_x_energy = log10(x_energy);
figure
hist(log_x_energy)
title('Histogram of signal log energy per frames')
