function[] = bad_channels(sFiles,badchan)

badchannels=badchan(1);
for i = 2:length(badchan)
    badchannels = strcat(badchannels,', ',badchan(i));

% Start a new report
bst_report('Start', sFiles);

% Process: Set bad channels
sFiles = bst_process('CallProcess', 'process_channel_setbad', sFiles, [], ...
    'sensortypes', badchannels);

% Save and display report
ReportFile = bst_report('Save', sFiles);
bst_report('Open', ReportFile);

