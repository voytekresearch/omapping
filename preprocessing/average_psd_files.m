function[]= average_psd_files(sFiles)


% Start a new report
bst_report('Start', sFiles);

% Process: Average: Everything
sFiles = bst_process('CallProcess', 'process_average', sFiles, [], ...
    'avgtype', 1, ...  % Everything
    'avg_func', 1, ...  % Arithmetic average:  mean(x)
    'weighted', 0, ...
    'matchrows', 1);

% Save and display report
ReportFile = bst_report('Save', sFiles);
bst_report('Open', ReportFile);

end