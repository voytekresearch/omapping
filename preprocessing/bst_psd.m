function [TF, FreqVector, Nwin, Messages] = bst_psd( F, sfreq, WinLength, WinOverlap, BadSegments )
% BST_PSD: Compute the PSD of a set of signals using Welch method

% @=============================================================================
% This function is part of the Brainstorm software:
% http://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2016 University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Francois Tadel, 2012-2014



% Parse inputs
if (nargin < 3) || isempty(WinLength) || (WinLength == 0)
    WinLength = size(F,2) ./ sfreq;
end
if (nargin < 4) || isempty(WinOverlap)
    WinOverlap = 50;
end
if (nargin < 5) || isempty(BadSegments)
    BadSegments = [];
end
Messages = '';
% Get sampling frequency
nTime = size(F,2);

% Initialize returned values
TF = [];

FreqVector = [];
Nwin = [];

% ===== WINDOWING =====
Lwin  = round(WinLength * sfreq);
Loverlap = round(Lwin * WinOverlap / 100);
% If window is too small
if (Lwin < 50)
    Messages = ['Time window is too small, please increase it and run the process again.' 10];
    return;
% If window is bigger than the data
elseif (Lwin > nTime)
    Lwin = size(F,2);
    Loverlap = 0;
    Nwin = 1;
    Messages = ['Time window is too large, using the entire recordings to estimate the spectrum.' 10];
% Else: there is at least one full time window
else
    Nwin = floor((nTime - Loverlap) ./ (Lwin - Loverlap));
    Messages = [Messages, sprintf('Using %d windows of %d samples each', Nwin, Lwin)];
end
% Next power of 2 from length of signal
NFFT = 2^nextpow2(Lwin);
% Positive frequency bins spanned by FFT
FreqVector = sfreq / 2 * linspace(0,1,NFFT/2+1);


% ===== CALCULATE FFT FOR EACH WINDOW =====
Nbad = 0;
i=1;
TF1=[];
for iWin = 1:Nwin
    % Build indices
    iTimes = (1:Lwin) + (iWin-1)*(Lwin - Loverlap);
    % Check if this segment is outside of ALL the bad segments (either entirely before or entirely after)
    if ~isempty(BadSegments) && (~all((iTimes(end) < BadSegments(1,:)) | (iTimes(1) > BadSegments(2,:))))
        disp(sprintf('BST> Skipping window #%d because it contains a bad segment.', iWin));
        Nbad = Nbad + 1;
        continue;
    end
    % Select indices
    Fwin = F(:,iTimes);
    % Remove mean of the signal
    Fwin = bst_bsxfun(@minus, Fwin, mean(Fwin,2));
    % Apply a hamming window to signal
    Fwin = bst_bsxfun(@times, Fwin, bst_window('hamming', Lwin)');
    % Compute FFT
    Ffft = fft(Fwin, NFFT, 2);
    % Keep only first half
    % (x2 to recover full power from negative frequencies)
    TFwin = 2 * Ffft(:,1:NFFT/2+1) ./ length(iTimes);
    z = TFwin;

    % Permute dimensions: time and frequency
    TFwin = permute(TFwin, [1 3 2]);
    % Convert to power
    TFwin = process_tf_measure('Compute', TFwin, 'none', 'power');
    % Add PSD of the window to the average
    TF1(:,:,:,i) = TFwin;
    i=i+1;
%     if isempty(TF)
%         TF = TFwin ./ Nwin;
%     else
%         TF = TF + TFwin ./ Nwin;
%     end 
end

[r c t f] = size(TF1);

for k = 1:t
    for j = 1:r
        TF(j,1,k)=median(TF1(j,1,k,1:f));
    end
end

    
% % Correct the dividing factor if there are bad segments
% if (Nbad > 0)
%     TF = TF .* (Nwin ./ (Nwin - Nbad));
%     Nwin = Nwin - Nbad;
% end





