function [mLogMelSpectrogram, vMelFrequencies, v_time] = ComputeLogMelSpectrogram(vSignal, fs, frame_length, frame_shift)
	[m_stft, v_freq, v_time] = my_stft(vSignal, fs, frame_length, frame_shift, hann(frame_length*fs, 'periodic'));
	mSTFTAmplitudeSquared = abs(m_stft).^2;
	load('MelFilterbank.mat');
	mMelSpectrogram = mMelFilterbank * mSTFTAmplitudeSquared;
	mLogMelSpectrogram = log10(mMelSpectrogram);
endfunction
