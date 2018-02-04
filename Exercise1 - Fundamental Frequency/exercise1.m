close all

[M_1,fs_1] = audioread ("AudioFiles/speech1.wav");
[M_2,fs_2] = audioread ("AudioFiles/speech2.wav");


[m_frames_1, v_time_frame_1] = my_windowing(M_1, fs_1, 0.032, 0.016);
[m_frames_2, v_time_frame_2] = my_windowing(M_2, fs_2, 0.032, 0.016);


for i=1:rows(m_frames_1)
	flipped_frame =  fliplr(m_frames_1(i,:));
	acf_1(i,1:length(m_frames_1(i,:))*2-1) = conv(m_frames_1(i,:), flipped_frame);
endfor


for i=1:rows(m_frames_2)
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
endfor

zero_idx = length(m_frames_2(1,:))+1;
start_idx = zero_idx + cast(1.0/400 * fs_2, "int16");
end_idx = zero_idx + cast(1.0/80 * fs_2, "int16");
for i=1:length(acf_2(:,1))
	[max_val, max_index] = max(acf_2(i, start_idx:end_idx));
	true_max_idx = start_idx+max_index-1-zero_idx;
	f_period_len = double(true_max_idx)/double(fs_2);
	f_freqs_2(i) = 1/f_period_len;
endfor


% plot the estimated pitches

t1= linspace(0, length(M_1)/fs_1, length(M_1));
t2= linspace(0, length(M_2)/fs_2, length(M_2));
x1= v_time_frame_1/fs_1;
x2= v_time_frame_2/fs_2;

figure
subplot(2,1,1)       % add first plot in 2 x 1 grid
plot(t1,M_1)
xlabel('time (s)')
title('Speech Signal 1')

subplot(2,1,2)       % add second plot in 2 x 1 grid
plot(x1,f_freqs_1)
xlabel('time (s)')
ylabel('Frequency (Hz)')
title('Estimated fundamental frequency')
drawnow

figure
subplot(2,1,1)       % add first plot in 2 x 1 grid
plot(t2,M_2)
xlabel('time (s)')
title('Speech Signal 2')

subplot(2,1,2)       % add second plot in 2 x 1 grid
plot(x2,f_freqs_2)
xlabel('time (s)')
ylabel('Frequency (Hz)')     
title('Estimated fundamental frequency')

