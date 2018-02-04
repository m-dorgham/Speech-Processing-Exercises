close all

pkg load signal

frame_length = 0.032;
frame_shift = 0.010;

[SigA,fsA] = audioread ("AudioData/stop_normal_df.wav");
[mLogMelSpectrogramA, vMelFrequenciesA, vTimeFrameA] = ComputeLogMelSpectrogram(SigA, fsA, frame_length, frame_shift);

figure
imagesc(vTimeFrameA/fsA, 1:23, mLogMelSpectrogramA);
axis xy;
set(gca, 'YTick', 1:23);
% Set the correct frequencies by exploiting the previously employed integer
% numbering for the Mel bands .
set(gca, 'YTickLabel', round(vMelFrequenciesA(get(gca, 'YTick'))));
colorbar
title('Log Mel Spectrogram for the word stop (normal) of speaker df')
drawnow

%section 3: DTW
[SigB,fsB] = audioread ("AudioData/stop_slow_df.wav");
[mLogMelSpectrogramB, vMelFrequenciesB, vTimeFrameB] = ComputeLogMelSpectrogram(SigB, fsB, frame_length, frame_shift);

[DTWCost, mDistance, vPath] = DynamicTimeWarp(mLogMelSpectrogramA, mLogMelSpectrogramB);

figure
imagesc(mDistance)
colormap('gray')
axis xy;
colorbar
title('Distance matrix between stop (normal and slow) of speaker df')
drawnow

%construct the warped Log-Mel-Spectrogram of the two signals
WmLogMelSpectrogramA(:,1) = mLogMelSpectrogramA(:,end);
WmLogMelSpectrogramB(:,1) = mLogMelSpectrogramB(:,end);
idxA=2;  idxB=2;
n=columns(mLogMelSpectrogramA);
m=columns(mLogMelSpectrogramB);
deltaTimeA = vTimeFrameA(2)-vTimeFrameA(1);
deltaTimeB = vTimeFrameB(2)-vTimeFrameB(1);

for i=1:length(vPath)
	if vPath(i) == 1
		WmLogMelSpectrogramA(:,idxA) = mLogMelSpectrogramA(:,n);
		idxA++;
		WmLogMelSpectrogramB(:,idxB) = mLogMelSpectrogramB(:,m-1);
		idxB++;  m--;
	elseif vPath(i) == 3
		WmLogMelSpectrogramA(:,idxA) = mLogMelSpectrogramA(:,n-1);
		idxA++;  n--;
		WmLogMelSpectrogramB(:,idxB) = mLogMelSpectrogramB(:,m);
		idxB++;  
	elseif vPath(i) == 2
		WmLogMelSpectrogramA(:,idxA) = mLogMelSpectrogramA(:,n-1);
		idxA++;  n--;
		WmLogMelSpectrogramB(:,idxB) = mLogMelSpectrogramB(:,m-1);
		idxB++;  m--;	
	end
end
%now reverse the warped spectrograms since they were constructed in a reversed order
WmLogMelSpectrogramA = fliplr(WmLogMelSpectrogramA);
WmLogMelSpectrogramB = fliplr(WmLogMelSpectrogramB);
WTimeFrameA = vTimeFrameA:deltaTimeA:columns(WmLogMelSpectrogramA)*deltaTimeA;
WTimeFrameB = vTimeFrameB:deltaTimeB:columns(WmLogMelSpectrogramB)*deltaTimeB;

figure
set(gcf, 'Position', [100, 200, 1100, 400])
subplot(1,2,1);
imagesc(WTimeFrameA/fsA, 1:23, WmLogMelSpectrogramA);
set(gca, 'CLim', [-8 4])
axis xy;
set(gca, 'YTickLabel', round(vMelFrequenciesA(get(gca, 'YTick'))));
colorbar
xlabel('time (s)')
ylabel('Mel band')
title('First signal warped Log-Mel-spectrogram')
subplot(1,2,2);
imagesc(WTimeFrameB/fsB, 1:23, WmLogMelSpectrogramB);
set(gca, 'CLim', [-8 4])
axis xy;
set(gca, 'YTickLabel', round(vMelFrequenciesB(get(gca, 'YTick'))));
colorbar
xlabel('time (s)')
ylabel('Mel band')
title('Second signal warped Log-Mel-spectrogram')
drawnow

%section 4: training on single speaker

AudioPath = 'AudioData';
stDfNormalFiles = dir(fullfile(AudioPath, '*_normal_df.wav'));
[vSignal_1, fs_1] = audioread(fullfile(AudioPath, stDfNormalFiles(1).name));
[vSignal_2, fs_2] = audioread(fullfile(AudioPath, stDfNormalFiles(2).name));
[vSignal_3, fs_3] = audioread(fullfile(AudioPath, stDfNormalFiles(3).name));
[vSignal_4, fs_4] = audioread(fullfile(AudioPath, stDfNormalFiles(4).name));
[vSignal_5, fs_5] = audioread(fullfile(AudioPath, stDfNormalFiles(5).name));


[mLogMelSpectrogram_1, vMelFrequencies_1, vTimeFrame_1] = ComputeLogMelSpectrogram(vSignal_1, fs_1, frame_length, frame_shift);
[mLogMelSpectrogram_2, vMelFrequencies_2, vTimeFrame_2] = ComputeLogMelSpectrogram(vSignal_2, fs_2, frame_length, frame_shift);
[mLogMelSpectrogram_3, vMelFrequencies_3, vTimeFrame_3] = ComputeLogMelSpectrogram(vSignal_3, fs_3, frame_length, frame_shift);
[mLogMelSpectrogram_4, vMelFrequencies_4, vTimeFrame_4] = ComputeLogMelSpectrogram(vSignal_4, fs_4, frame_length, frame_shift);
[mLogMelSpectrogram_5, vMelFrequencies_5, vTimeFrame_5] = ComputeLogMelSpectrogram(vSignal_5, fs_5, frame_length, frame_shift);

DTWCosts = zeros(5,5);

[DTWCost, mDistance, vPath] = DynamicTimeWarp(mLogMelSpectrogram_1, mLogMelSpectrogram_1);
DTWCosts(1,1) = DTWCost;
[DTWCost, mDistance, vPath] = DynamicTimeWarp(mLogMelSpectrogram_2, mLogMelSpectrogram_2);
DTWCosts(2,2) = DTWCost;
[DTWCost, mDistance, vPath] = DynamicTimeWarp(mLogMelSpectrogram_3, mLogMelSpectrogram_3);
DTWCosts(3,3) = DTWCost;
[DTWCost, mDistance, vPath] = DynamicTimeWarp(mLogMelSpectrogram_4, mLogMelSpectrogram_4);
DTWCosts(4,4) = DTWCost;
[DTWCost, mDistance, vPath] = DynamicTimeWarp(mLogMelSpectrogram_5, mLogMelSpectrogram_5);
DTWCosts(5,5) = DTWCost;
[DTWCost, mDistance, vPath] = DynamicTimeWarp(mLogMelSpectrogram_1, mLogMelSpectrogram_2);
DTWCosts(1,2) = DTWCost;	DTWCosts(2,1) = DTWCost;
[DTWCost, mDistance, vPath] = DynamicTimeWarp(mLogMelSpectrogram_1, mLogMelSpectrogram_3);
DTWCosts(1,3) = DTWCost;	DTWCosts(3,1) = DTWCost;
[DTWCost, mDistance, vPath] = DynamicTimeWarp(mLogMelSpectrogram_1, mLogMelSpectrogram_4);
DTWCosts(1,4) = DTWCost;	DTWCosts(4,1) = DTWCost;
[DTWCost, mDistance, vPath] = DynamicTimeWarp(mLogMelSpectrogram_1, mLogMelSpectrogram_5);
DTWCosts(1,5) = DTWCost;	DTWCosts(5,1) = DTWCost;
[DTWCost, mDistance, vPath] = DynamicTimeWarp(mLogMelSpectrogram_2, mLogMelSpectrogram_3);
DTWCosts(2,3) = DTWCost;	DTWCosts(3,2) = DTWCost;
[DTWCost, mDistance, vPath] = DynamicTimeWarp(mLogMelSpectrogram_2, mLogMelSpectrogram_4);
DTWCosts(2,4) = DTWCost;	DTWCosts(4,2) = DTWCost;
[DTWCost, mDistance, vPath] = DynamicTimeWarp(mLogMelSpectrogram_2, mLogMelSpectrogram_5);
DTWCosts(2,5) = DTWCost;	DTWCosts(5,2) = DTWCost;
[DTWCost, mDistance, vPath] = DynamicTimeWarp(mLogMelSpectrogram_3, mLogMelSpectrogram_4);
DTWCosts(3,4) = DTWCost;	DTWCosts(4,3) = DTWCost;
[DTWCost, mDistance, vPath] = DynamicTimeWarp(mLogMelSpectrogram_3, mLogMelSpectrogram_5);
DTWCosts(3,5) = DTWCost;	DTWCosts(5,3) = DTWCost;
[DTWCost, mDistance, vPath] = DynamicTimeWarp(mLogMelSpectrogram_4, mLogMelSpectrogram_5);
DTWCosts(4,5) = DTWCost;	DTWCosts(5,4) = DTWCost;

for i = 1:5
	cSplits = regexp(stDfNormalFiles(i).name(1:end-4), '_', 'split');
	tickLabels{i} = cSplits{1};
end

figure
imagesc(DTWCosts)
set(gca, 'XTickLabel', tickLabels)
set(gca, 'YTickLabel', tickLabels)
colormap('gray')
colorbar
title("DTW distance of different training words.")
drawnow


%section 5: recognition experiments

disp('Recognizing slow utterances of speaker df:')
stDfSlowFiles = dir(fullfile(AudioPath, '*_slow_df.wav'));
for i=1:length(stDfSlowFiles)
	RecognizedWord = RecognizeDTW(AudioPath, stDfNormalFiles, stDfSlowFiles(i));
	cSplits = regexp(stDfSlowFiles(i).name(1:end-4), '_', 'split');
	word = cSplits{1};
	str = ['Test Word: ', word, ',	Recognized Word: ', RecognizedWord];
	disp(str);
endfor

disp('*******************************************')
disp('Recognizing all utterances of speaker bj against df:')
stDfAllFiles = dir(fullfile(AudioPath, '*_df.wav'));
stBjAllFiles = dir(fullfile(AudioPath, '*_bj.wav'));
mistakes = 0;
for i=1:length(stBjAllFiles)
	RecognizedWord = RecognizeDTW(AudioPath, stDfAllFiles, stBjAllFiles(i));
	cSplits = regexp(stBjAllFiles(i).name(1:end-4), '_', 'split');
	word = cSplits{1};	speed = cSplits{2};
	str = ['Test Word: ', word, ' (', speed, '),	Recognized Word: ', RecognizedWord];
	disp(str)
	if ~strcmpi(word, RecognizedWord) 
		mistakes++;
	end
endfor
err_rate = mistakes/length(stBjAllFiles);
str = ['error rate for speaker bj (compared to df): ', num2str(err_rate)];
disp(str);

disp('*******************************************')
disp('Recognizing all utterances of speaker sw against df:')
stDfAllFiles = dir(fullfile(AudioPath, '*_df.wav'));
stSwAllFiles = dir(fullfile(AudioPath, '*_sw.wav'));
mistakes = 0;
for i=1:length(stSwAllFiles)
	RecognizedWord = RecognizeDTW(AudioPath, stDfAllFiles, stSwAllFiles(i));
	cSplits = regexp(stSwAllFiles(i).name(1:end-4), '_', 'split');
	word = cSplits{1};	speed = cSplits{2};
	str = ['Test Word: ', word, ' (', speed, '),	Recognized Word: ', RecognizedWord];
	disp(str)
	if ~strcmpi(word, RecognizedWord) 
		mistakes++;
	end
endfor
err_rate = mistakes/length(stSwAllFiles);
str = ['error rate for speaker sw (compared to df): ', num2str(err_rate)];
disp(str);

disp('*******************************************')
disp('Recognizing all utterances of speaker bj against sw:')
stSwAllFiles = dir(fullfile(AudioPath, '*_sw.wav'));
stBjAllFiles = dir(fullfile(AudioPath, '*_bj.wav'));
mistakes = 0;
for i=1:length(stBjAllFiles)
	RecognizedWord = RecognizeDTW(AudioPath, stSwAllFiles, stBjAllFiles(i));
	cSplits = regexp(stBjAllFiles(i).name(1:end-4), '_', 'split');
	word = cSplits{1};	speed = cSplits{2};
	str = ['Test Word: ', word, ' (', speed, '),	Recognized Word: ', RecognizedWord];
	disp(str)
	if ~strcmpi(word, RecognizedWord) 
		mistakes++;
	end
endfor
err_rate = mistakes/length(stBjAllFiles);
str = ['error rate for speaker bj (compared to sw): ', num2str(err_rate)];
disp(str);

disp('*******************************************')
disp('Recognizing all utterances of speaker sk against sw:')
stSwAllFiles = dir(fullfile(AudioPath, '*_df.wav'));
stSkAllFiles = dir(fullfile(AudioPath, '*_sw.wav'));
mistakes = 0;
for i=1:length(stSkAllFiles)
	RecognizedWord = RecognizeDTW(AudioPath, stSwAllFiles, stSkAllFiles(i));
	cSplits = regexp(stSkAllFiles(i).name(1:end-4), '_', 'split');
	word = cSplits{1};	speed = cSplits{2};
	str = ['Test Word: ', word, ' (', speed, '),	Recognized Word: ', RecognizedWord];
	disp(str)
	if ~strcmpi(word, RecognizedWord) 
		mistakes++;
	end
endfor
err_rate = mistakes/length(stSkAllFiles);
str = ['error rate for speaker sk (compared to sw): ', num2str(err_rate)];
disp(str);

