function[] = dc_offest_block(sFiles)


% Start a new report
bst_report('Start', sFiles);

% Process: Remove DC offset: [All file]
sFiles = bst_process('CallProcess', 'process_baseline', sFiles, [], ...
    'baseline', [], ...
    'sensortypes', 'MEG', ...
    'overwrite', 1);

% Save and display report
ReportFile = bst_report('Save', sFiles);
bst_report('Open', ReportFile);

end