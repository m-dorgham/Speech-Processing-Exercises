close all
pkg load signal

[S_w,fs_w] = audioread ("Audio/SpeechWhite.wav");
[S_b,fs_b] = audioread ("Audio/SpeechBabble.wav");

frame_length = 0.032;
frame_shift = 0.016;
theta=15;

[m_stft_w, v_freq_w, v_time_w] = my_stft(S_w, fs_w, frame_length, frame_shift, sqrt(hann(frame_length*fs_w, 'periodic')));
[m_stft_b, v_freq_b, v_time_b] = my_stft(S_b, fs_b, frame_length, frame_shift, sqrt(hann(frame_length*fs_b, 'periodic')));


%initial noise periodogram - use the first frame
noise_w_PSD(:,1) = abs(m_stft_w(:,1)).^2;

Q_w(:,1) = zeros(length(v_freq_w),1);
spp_w = zeros(length(v_freq_w), columns(m_stft_w));

S_w_PSD = abs(m_stft_w).^2;

for i=1:columns(S_w_PSD)
	spp_w(:,i) = (1+(1+theta)*exp(-(S_w_PSD(:,i)./noise_w_PSD(:,i))*(theta/(1+theta)))).^-1;
	Q_w(:,i+1) = 0.9*Q_w(:,i) + 0.1*spp_w(:,i);
	for j=1:rows(Q_w)
		if(Q_w(j,i+1)>0.99)
			spp_w(j,i) = min(0.99, spp_w(j,i));
		endif
	endfor
	N_w_periodogram = spp_w(:,i).*noise_w_PSD(:,i) + (1-spp_w(:,i)).*S_w_PSD(:,i);
	noise_w_PSD(:,i+1) = 0.8*noise_w_PSD(:,i) + 0.2*N_w_periodogram;
endfor


figure
imagesc(v_time_w/fs_w, v_freq_w, spp_w)
axis xy
colorbar
xlabel('time (s)')
ylabel('Frequency (Hz)')
title('Speech presence probability spectroram')
figure
imagesc(v_time_w/fs_w, v_freq_w, 10*log10(noise_w_PSD(:,2:end)))
axis xy
colorbar
xlabel('time (s)')
ylabel('Frequency (Hz)')
title('Noise PSD spectrogram')
figure
imagesc(v_time_w/fs_w, v_freq_w, 10*log10(S_w_PSD))
axis xy
colorbar
xlabel('time (s)')
ylabel('Frequency (Hz)')
title('Speech signal spectrogram')

%repeat for speech babble signal
noise_b_PSD(:,1) = abs(m_stft_b(:,1)).^2;

Q_b(:,1) = zeros(length(v_freq_b),1);
spp_b = zeros(length(v_freq_b), columns(m_stft_b));

S_b_PSD = abs(m_stft_b).^2;

for i=1:columns(S_b_PSD)
	spp_b(:,i) = (1+(1+theta)*exp(-(S_b_PSD(:,i)./noise_b_PSD(:,i))*(theta/(1+theta)))).^-1;
	Q_b(:,i+1) = 0.9*Q_b(:,i) + 0.1*spp_b(:,i);
	for j=1:rows(Q_b)
		if(Q_b(j,i+1)>0.99)
			spp_b(j,i) = min(0.99, spp_b(j,i));
		endif
	endfor
	N_b_periodogram = spp_b(:,i).*noise_b_PSD(:,i) + (1-spp_b(:,i)).*S_b_PSD(:,i);
	noise_b_PSD(:,i+1) = 0.8*noise_b_PSD(:,i) + 0.2*N_b_periodogram;
endfor

figure
imagesc(v_time_b/fs_b, v_freq_b, spp_b)
axis xy
colorbar
xlabel('time (s)')
ylabel('Frequency (Hz)')
title('Speech presence probability spectroram')
figure
imagesc(v_time_b/fs_b, v_freq_b, 10*log10(noise_b_PSD(:,2:end)))
axis xy
colorbar
xlabel('time (s)')
ylabel('Frequency (Hz)')
title('Noise PSD spectrogram')
figure
imagesc(v_time_b/fs_b, v_freq_b, 10*log10(S_b_PSD))
axis xy
colorbar
xlabel('time (s)')
ylabel('Frequency (Hz)')
title('Speech signal spectrogram')


%A priori SNR estimation and Wiener Filtering
alpha = 0.5;
G_min = 0;

%repeat the previous process to compute the noise PSD and then compute the enhanced speech
noise_w_PSD(:,1) = abs(m_stft_w(:,1)).^2;

Q_w(:,1) = zeros(length(v_freq_w),1);
spp_w = zeros(length(v_freq_w), columns(m_stft_w));
matEnhancedSpeech_w(:,1) = zeros(length(v_freq_w),1);
zeta_w = zeros(length(v_freq_w), columns(m_stft_w));
G_w = zeros(length(v_freq_w), columns(m_stft_w));

S_w_PSD = abs(m_stft_w).^2;

for i=1:columns(S_w_PSD)
	spp_w(:,i) = (1+(1+theta)*exp(-(S_w_PSD(:,i)./noise_w_PSD(:,i))*(theta/(1+theta)))).^-1;
	Q_w(:,i+1) = 0.9*Q_w(:,i) + 0.1*spp_w(:,i);
	for j=1:rows(Q_w)
		if(Q_w(j,i+1)>0.99)
			spp_w(j,i) = min(0.99, spp_w(j,i));
		endif
	endfor
	N_w_periodogram = spp_w(:,i).*noise_w_PSD(:,i) + (1-spp_w(:,i)).*S_w_PSD(:,i);
	noise_w_PSD(:,i+1) = 0.8*noise_w_PSD(:,i) + 0.2*N_w_periodogram;

	zeta_w(:,i) = alpha*matEnhancedSpeech_w(:,i)./noise_w_PSD(:,i) + (1-alpha)*max((S_w_PSD(:,i)./noise_w_PSD(:,i+1))-1,0);
	G_w(:,i) = max(zeta_w(:,i)./(1+zeta_w(:,i)),G_min);
	
	matEnhancedSpeech_w(:,i+1) = G_w(:,i).*m_stft_w(:,i);
endfor


figure
set(gcf, 'Position', [100, 200, 1100, 400])
subplot(1,2,1);
imagesc(v_time_w/fs_w, v_freq_w, 10*log10(max(abs(m_stft_w).^2, 10^(-15))))
set(gca, 'CLim', [-100 40])
axis xy
colorbar
xlabel('time (s)')
ylabel('Frequency (Hz)')
title('Noisy speech signal spectrogram')
subplot(1,2,2);
imagesc(v_time_w/fs_w, v_freq_w, 10*log10(max(abs(matEnhancedSpeech_w(:,2:end)).^2, 10^(-15))))
set(gca, 'CLim', [-100 40])
axis xy
colorbar
xlabel('time (s)')
ylabel('Frequency (Hz)')
title('Enhanced speech signal spectrogram')

%repeat for babble noisy signal
noise_b_PSD(:,1) = abs(m_stft_b(:,1)).^2;

Q_b(:,1) = zeros(length(v_freq_b),1);
spp_b = zeros(length(v_freq_b), columns(m_stft_b));
matEnhancedSpeech_b(:,1) = zeros(length(v_freq_b),1);
zeta_b = zeros(length(v_freq_b), columns(m_stft_b));
G_b = zeros(length(v_freq_b), columns(m_stft_b));

S_b_PSD = abs(m_stft_b).^2;

for i=1:columns(S_b_PSD)
	spp_b(:,i) = (1+(1+theta)*exp(-(S_b_PSD(:,i)./noise_b_PSD(:,i))*(theta/(1+theta)))).^-1;
	Q_b(:,i+1) = 0.9*Q_b(:,i) + 0.1*spp_b(:,i);
	for j=1:rows(Q_b)
		if(Q_b(j,i+1)>0.99)
			spp_b(j,i) = min(0.99, spp_b(j,i));
		endif
	endfor
	N_b_periodogram = spp_b(:,i).*noise_b_PSD(:,i) + (1-spp_b(:,i)).*S_b_PSD(:,i);
	noise_b_PSD(:,i+1) = 0.8*noise_b_PSD(:,i) + 0.2*N_b_periodogram;

	zeta_b(:,i) = alpha*matEnhancedSpeech_b(:,i)./noise_b_PSD(:,i) + (1-alpha)*max((S_b_PSD(:,i)./noise_b_PSD(:,i+1))-1,0);
	G_b(:,i) = max(zeta_b(:,i)./(1+zeta_b(:,i)),G_min);
	
	matEnhancedSpeech_b(:,i+1) = G_b(:,i).*m_stft_b(:,i);
endfor


figure
set(gcf, 'Position', [100, 200, 1100, 400])
subplot(1,2,1);
imagesc(v_time_b/fs_b, v_freq_b, 10*log10(max(abs(m_stft_b).^2, 10^(-15))))
set(gca, 'CLim', [-100 40])
axis xy
colorbar
xlabel('time (s)')
ylabel('Frequency (Hz)')
title('Noisy speech signal spectrogram')
subplot(1,2,2);
imagesc(v_time_b/fs_b, v_freq_b, 10*log10(max(abs(matEnhancedSpeech_b(:,2:end)).^2, 10^(-15))))
set(gca, 'CLim', [-100 40])
axis xy
colorbar
xlabel('time (s)')
ylabel('Frequency (Hz)')
title('Enhanced speech signal spectrogram')


synthesizedEnhancedSignal_w = my_inverse_stft(matEnhancedSpeech_w(:,2:end), fs_w, frame_length, frame_shift, sqrt(hann(frame_length*fs_w, 'periodic')));
synthesizedEnhancedSignal_b = my_inverse_stft(matEnhancedSpeech_b(:,2:end), fs_b, frame_length, frame_shift, sqrt(hann(frame_length*fs_b, 'periodic')));

soundsc(S_w, fs_w)
soundsc(synthesizedEnhancedSignal_w, fs_w)
soundsc(S_b, fs_b)
soundsc(synthesizedEnhancedSignal_b, fs_b)
