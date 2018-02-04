% Computes inverse STFT from a STFT matrix
% v_signal = my_inverse_stft(m_stft, fs, frame_length, frame_shift, v_synthesis_window)
%
%   Inputs:
%       m_stft              -       STFT matrix
%       fs                  -       sampling rate in Hz
%       frame_length        -       frame length in seconds that was used for separating the signal 
%                                   into blocks
%       frame_shift         -       frame shift in seconds used for the frames in m_stft
%       v_synthesis_window  -       vector containing a synthesis window function (needs to have 
%                                   same length as the frame length in samples)
%
%   Outputs:
%       v_signal            -       signal synthesized from STFT

function v_signal = my_inverse_stft(m_stft, fs, frame_length, frame_shift, v_synthesis_window)

% File   :  my_inverse_stft.m
% Author :  Robert Rehr <r.rehr AT uni-oldenburg.de>
% Date   :  02.11.2017
%
% Updates:

% get number of frames
num_frames = size(m_stft, 2);

% compute fft length
fft_length = round(fs * frame_length);

% compute frame_shift in samples
frame_shift = round(frame_shift * fs);

% restore complex conjugate spectrum
m_stft(ceil(fft_length / 2) + 1:fft_length, :) = conj(m_stft(end:-1:2, :));

% compute inverse stft and apply synthesis window
m_frames = bsxfun(@times, real(ifft(m_stft)), v_synthesis_window(:));

% compute indeces for each frame
m_idx = bsxfun(@plus, (1:fft_length).', (0:num_frames - 1) * frame_shift);

% reconstruct time domain signal
v_signal = accumarray(m_idx(:), m_frames(:));

% End of my_inverse_stft.m
