function[] = dc_offest_block(sFiles)


% Start a new report
bst_report('Start', sFiles);

% Process: DC offset correction: [All file]
sFiles = bst_process('CallProcess', 'process_baseline', sFiles, [], ...
    'baseline',    [], ...
    'sensortypes', 'MEG', ...
    'method',      'bl', ...  % DC offset correction:    x_std = x - &mu;
    'overwrite',   1);

% Save and display report
ReportFile = bst_report('Save', sFiles);
bst_report('Open', ReportFile);

end