function [sigOut, filterStateOut] = filterAdaptively(MACoeff, ARCoeff, sigIn, filterStateIn)
% Allows for segment-wise filtering of a signal with changing filters.
% Input:  MACoeff        - Moving average filter coefficients (vector)
%         ARCoeff        - Autoregressive filter coefficients (e.g. LPCs) (vector)
%         sigIN          - Input signal segment (vector)
%         filterStateIn  - 'State' of the filter before filtering
% 
% Output: sigOut         - Output signal segment (filtered version of sigIn)
%         filterStateOut - 'State' of the filter after filtering
%
% Usage:  Example for LPC filtering:
%         Call [segmentOut,filtState] = filterAdaptively(1, LPCs,segmentIn,filtState);
%         for every signal segment, using the corresponding (time varying)
%         LPCs for this frame. 'filterAdaptively' will ensure a correct
%         initialization of the time varying filter for each segment.
%         For the first segment, do not use 'filtState' as an input,
%         'filterAdaptively' will then initialize and return the first
%         filter state.

    sigIn = sigIn(:);
    shift = length(sigIn);

    if nargin < 4 % Create new filter state if none is provided
        filterStateIn = [];
    end

    if isempty(filterStateIn) % Initialize new filter state
        filterStateIn.prevIn  = zeros(length(MACoeff)-1,1);
        filterStateIn.prevOut = zeros(length(ARCoeff)-1,1);
    end

    filterStateInternal    = filtic(MACoeff, ARCoeff, filterStateIn.prevOut, filterStateIn.prevIn); % Set initial State
    sigOut                 = filter(MACoeff, ARCoeff, sigIn, filterStateInternal); % Apply filter

    filterStateOut.prevIn  = [flipud(sigIn);filterStateIn.prevIn(1:end-shift)];   % Return current filter state
    filterStateOut.prevOut = [flipud(sigOut);filterStateIn.prevOut(1:end-shift)]; % Return current filter state

end