function[] = resample(sFiles)

% Start a new report
bst_report('Start', sFiles);

% Process: Resample: 1000Hz
sFiles = bst_process('CallProcess', 'process_resample', sFiles, [], ...
    'freq',      1000, ...
    'overwrite', 1);

% Save and display report
ReportFile = bst_report('Save', sFiles);
bst_report('Open', ReportFile);

