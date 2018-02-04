close all
pkg load signal

[sig,fs] = audioread ("Audio/speech1.wav");

%plot(sig)

frame_len_samples = 0.032*fs;

voiced_frame = sig(5400:5400+frame_len_samples-1); %the starting index is obtained by ploting the signal
unvoiced_frame = sig(16000:16000+frame_len_samples-1);

windowed_voiced_frame = voiced_frame.*hann(frame_len_samples, 'periodic');
windowed_unvoiced_frame = unvoiced_frame.*hann(frame_len_samples, 'periodic');

lpc_order=12;

%compute the autocorrelation vector
phi_v = xcorr(windowed_voiced_frame);
phi_u = xcorr(windowed_unvoiced_frame);

phi_v = phi_v(round(length(phi_v)/2):round(length(phi_v)/2)+lpc_order);
phi_u = phi_u(round(length(phi_u)/2):round(length(phi_u)/2)+lpc_order);

%compute the toeplitz matrix
R_v = toeplitz(phi_v(1:end-1));
R_u = toeplitz(phi_u(1:end-1));

phi_v = phi_v(2:end);
phi_u = phi_u(2:end);

%compute LP coefficients
a_v = -(R_v\phi_v);
a_u = -(R_u\phi_u);

%frequency response of the estimated vocal tract filter
[H_v, freq_ax1] = freqz(1,[1;a_v],frame_len_samples,'whole',fs);
[H_u, freq_ax2] = freqz(1,[1;a_u],frame_len_samples,'whole',fs);

%compute the DFT of the windowed segments
S_v = fft(windowed_voiced_frame);
S_u = fft(windowed_unvoiced_frame);
sigma_sv = sqrt(norm((S_v))^2/length(S_v));
sigma_su = sqrt(norm((S_u))^2/length(S_u));
sigma_hv = sqrt(norm((H_v))^2/length(H_v));
sigma_hu = sqrt(norm((H_u))^2/length(H_u));
%normalize H to align its plot on S plot
H_v = H_v.*sigma_sv/sigma_hv;
H_u = H_u.*sigma_su/sigma_hu;

figure
plot(freq_ax1, 20*log10(abs(H_v)))
hold on
plot(freq_ax1, 20*log10(abs(S_v)))
hold off
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
title('frequency response of the estimated vocal tract filter (voiced segment)')
legend('Estimated response', 'Original response', 'Location','southeast')
axis tight
drawnow

figure
plot(freq_ax2, 20*log10(abs(H_u)))
hold on
plot(freq_ax2, 20*log10(abs(S_u)))
hold off
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
title('frequency response of the estimated vocal tract filter (unvoiced segment)')
legend('Estimated response', 'Original response', 'Location','southeast')
axis tight
drawnow

e_v = filter([1; a_v],1,voiced_frame);
e_u = filter([1; a_u],1,unvoiced_frame);

figure
plot(1:length(voiced_frame), voiced_frame, 1:length(e_v), e_v, '--')
legend('Original signal','Residual Signal')
title('The residual signal (voiced segment)')
drawnow

figure
plot(1:length(unvoiced_frame), unvoiced_frame, 1:length(e_u), e_u, '--')
legend('Original signal','Residual Signal')
title('The residual signal (unvoiced segment)')

