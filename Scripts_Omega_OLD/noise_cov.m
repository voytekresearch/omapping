function[] = noise_cov(sFiles)


% Start a new report
bst_report('Start', sFiles);

% Process: Compute noise covariance
sFiles = bst_process('CallProcess', 'process_noisecov', sFiles, [], ...
    'baseline', [], ...
    'sensortypes', 'MEG, EEG, SEEG, ECOG', ...
    'target', 1, ...  % Noise covariance
    'dcoffset', 1, ...  % Block by block, to avoid effects of slow shifts in data
    'method', 1, ...  % Full noise covariance matrix
    'copycond', 1, ...
    'copysubj', 0, ...
    'replacefile', 1);  % Replace

% Save and display report
ReportFile = bst_report('Save', sFiles);
bst_report('Open', ReportFile);

end