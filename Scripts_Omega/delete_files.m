function[] = delete_files(sFiles)

% Start a new report
bst_report('Start', sFiles);

% Process: Delete selected files
sFiles = bst_process('CallProcess', 'process_delete', sFiles, [], ...
    'target', 1);  % Delete selected files

% Save and display report
ReportFile = bst_report('Save', sFiles);
bst_report('Open', ReportFile);

end