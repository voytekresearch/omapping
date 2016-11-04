function[]=link_raw(SubjectNames, RawFiles)

% Input files
sFiles = [];

% Start a new report
bst_report('Start', sFiles);

% Process: Create link to raw file
for i = 1:length(RawFiles)
    sFiles = bst_process('CallProcess', 'process_import_data_raw', sFiles, [], ...
        'subjectname', SubjectNames{1}, ...
        'datafile', {RawFiles{i}, 'CTF'}, ...
        'channelreplace', 0, ...
        'channelalign', 1, ...
        'evtmode', 'value');
end

% Save and display report
ReportFile = bst_report('Save', sFiles);
bst_report('Open', ReportFile);

end