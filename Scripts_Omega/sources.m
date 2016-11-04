function[] = sources(sFiles)

% Start a new report
bst_report('Start', sFiles);

% Process: Compute sources
sFiles = bst_process('CallProcess', 'process_inverse', sFiles, [], ...
    'Comment',     '', ...
    'method',      1, ...  % Minimum norm estimates (wMNE)
    'wmne',        struct(...
         'NoiseCov', [], ...
         'InverseMethod', 'wmne', ...
         'ChannelTypes', {{}}, ...
         'SNR', 3, ...
         'diagnoise', 0, ...
         'SourceOrient', {{'fixed'}}, ...
         'loose', 0.2, ...
         'depth', 1, ...
         'weightexp', 0.5, ...
         'weightlimit', 10, ...
         'regnoise', 1, ...
         'magreg', 0.1, ...
         'gradreg', 0.1, ...
         'eegreg', 0.1, ...
         'ecogreg', 0.1, ...
         'seegreg', 0.1, ...
         'fMRI', [], ...
         'fMRIthresh', [], ...
         'fMRIoff', 0.1, ...
         'pca', 1), ...
    'sensortypes', 'MEG', ...
    'output',      1);  % Kernel only: shared

% Save and display report
ReportFile = bst_report('Save', sFiles);
bst_report('Open', ReportFile);

