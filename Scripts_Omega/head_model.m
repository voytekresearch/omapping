function[] = head_model(sFiles)

% Start a new report
bst_report('Start', sFiles);

% Process: Compute head model
sFiles = bst_process('CallProcess', 'process_headmodel', sFiles, [], ...
    'Comment',     '', ...
    'sourcespace', 1, ...  % Cortex surface
    'volumegrid',  [], ...
    'meg',         3, ...  % Overlapping spheres
    'eeg',         1, ...  % 
    'ecog',        1, ...  % 
    'seeg',        1, ...  % 
    'openmeeg',    struct(...
         'BemFiles', {{}}, ...
         'BemNames', {{'Scalp', 'Skull', 'Brain'}}, ...
         'BemCond', [1, 0.0125, 1], ...
         'BemSelect', [1, 1, 1], ...
         'isAdjoint', 0, ...
         'isAdaptative', 1, ...
         'isSplit', 0, ...
         'SplitLength', 4000));

% Save and display report
ReportFile = bst_report('Save', sFiles);
bst_report('Open', ReportFile);

