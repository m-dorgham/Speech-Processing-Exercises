function RecognizedWord = RecognizeDTW ( AudioPath , stTrainFiles , stTestFile )
	[vSignal_t, fs_t] = audioread(fullfile(AudioPath, stTestFile.name));
	[mLogMelSpectrogram_t, vMelFrequencies_t, vTimeFrame_t] = ComputeLogMelSpectrogram(vSignal_t, fs_t, 0.032, 0.010);
	leastCost = inf;
	nearst_word_idx = -1;

	for i=1:length(stTrainFiles)
		[vSignal_i, fs_i] = audioread(fullfile(AudioPath, stTrainFiles(i).name));
		[mLogMelSpectrogram_i, vMelFrequencies_i, vTimeFrame_i] = ComputeLogMelSpectrogram(vSignal_i, fs_i, 0.032, 0.010);
		[DTWCost, mDistance, vPath] = DynamicTimeWarp(mLogMelSpectrogram_t, mLogMelSpectrogram_i);
		if DTWCost < leastCost
			leastCost = DTWCost;
			nearst_word_idx = i;
		endif
	endfor

	cSplits = regexp(stTrainFiles(nearst_word_idx).name(1:end-4), '_', 'split');
	RecognizedWord = cSplits{1};
endfunction
