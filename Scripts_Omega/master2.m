% This script should be run after master1 and bad channels have been
% marked.

% Automatically detects artifacts like heartbeats, eyeblinks etc.
pipelinecount=3;
sFiles=inputfiles(pipelinecount);
artifact_detection(sFiles);
pipelinecount=pipelinecount+1;

% Manually check for artifacts and then run master3