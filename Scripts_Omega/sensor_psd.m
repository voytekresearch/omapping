function[]=sensor_psd(sFiles)

% Start a new report
bst_report('Start', sFiles);

% Process: Power spectrum density (Welch)
sFiles = bst_process('CallProcess', 'process_psd', sFiles, [], ...
    'timewindow', [], ...
    'win_length', 2, ...
    'win_overlap', 50, ...
    'sensortypes', 'MEG, EEG', ...
    'edit', []);

% Save and display report
ReportFile = bst_report('Save', sFiles);
bst_report('Open', ReportFile);

end