function[] = notch_filter(sFiles)


% Start a new report
bst_report('Start', sFiles);

% Process: Notch filter: 60Hz 120Hz 180Hz 240Hz 300Hz 360Hz 420Hz 480Hz
sFiles = bst_process('CallProcess', 'process_notch', sFiles, [], ...
    'freqlist', [60, 120, 180, 240, 300, 360, 420, 480], ...
    'sensortypes', 'MEG, EEG', ...
    'read_all', 1);

% Save and display report
ReportFile = bst_report('Save', sFiles);
bst_report('Open', ReportFile);

