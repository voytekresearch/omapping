function[] = delete_folder(sFiles)


% Start a new report
bst_report('Start', sFiles);

% Process: Delete folders
sFiles = bst_process('CallProcess', 'process_delete', sFiles, [], ...
    'target', 2);  % Delete folders

% Save and display report
ReportFile = bst_report('Save', sFiles);
bst_report('Open', ReportFile);

end

