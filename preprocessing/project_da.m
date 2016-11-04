function[] = project_da(sFiles)

% Start a new report
bst_report('Start', sFiles);

% Process: Project on default anatomy
sFiles = bst_process('CallProcess', 'process_project_sources', sFiles, []);

% Save and display report
ReportFile = bst_report('Save', sFiles);
bst_report('Open', ReportFile);

end