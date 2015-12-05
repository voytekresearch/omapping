% Script generated by Brainstorm (24-Nov-2015)

% Input files
sFiles = {...
    'Subject01/@rawsubj002_spontaneous_20111102_01_AUX/data_0raw_subj002_spontaneous_20111102_01_AUX.mat', ...
    'Subject01/@rawsubj002_spontaneous_20111102_02_AUX/data_0raw_subj002_spontaneous_20111102_02_AUX.mat'};

% Start a new report
bst_report('Start', sFiles);

% Process: Power spectrum density (Welch)
sFiles = bst_process('CallProcess', 'process_psd', sFiles, [], ...
    'timewindow', [], ...
    'win_length', 1, ...
    'win_overlap', 50, ...
    'sensortypes', 'MEG, EEG', ...
    'edit', []);

% Save and display report
ReportFile = bst_report('Save', sFiles);
bst_report('Open', ReportFile);

