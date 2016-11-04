% This script should be run after anatomies are imported in Brainstorm and
% links created to RAW files

% Denotes step in pipeline
pipelinecount=1;

% Notch Filter (60Hz and harmonics) MEG resting state recording and empty room noise
sFiles=inputfiles(pipelinecount);
notch_filter(sFiles);
pipelinecount=pipelinecount+1;

% Calculate PSDs to find channel outliers
sFiles=inputfiles(pipelinecount);
psd_chanout(sFiles);
pipelinecount=pipelinecount+1;

% Manually delete bad channels and then run master2




