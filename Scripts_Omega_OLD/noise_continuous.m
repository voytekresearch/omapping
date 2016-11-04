function[] = noise_continuous(sFiles)


% Start a new report
bst_report('Start', sFiles);

% Process: Convert to continuous (CTF): Continuous
sFiles = bst_process('CallProcess', 'process_ctf_convert', sFiles, [], ...
    'rectype', 2);  % Continuous

% Save and display report
ReportFile = bst_report('Save', sFiles);
bst_report('Open', ReportFile);

