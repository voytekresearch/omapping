% Load PSD data and set to visualize
% TODO: 
%   1) Currently assumes default oscillation bands. 
%       Should update to use whatever bands are found in load file.
%   2) Figure out channel data plotting.

%% Set Up

% Load a list of variable names for a Bst file
load('var_names.mat')

% Set paths
datapath = '/Users/tom/Documents/Research/1-Projects/MEGmapping/2-Data/ExternalData/MEG/Viz/';
bst_path = '/Users/tom/brainstorm_db/om_viz/data/om_viz_subj/Visualizing-MEG-FOOF/';

%% Load FOOOF Data - Oscillation Frequencies

subj = 559176;
filename = [num2str(subj), '_Foof_Viz.mat'];

% Load FOOOF Data
load([datapath, filename]);
%clear filename datapath

%
fooof_dat = {slopes, theta, alpha, beta, lowgamma};
fooof_labels = {['Slopes_', num2str(subj)], ['Thetas_', num2str(subj)], ...
    ['Alphas_', num2str(subj)], ['Betas_', num2str(subj)]};

disp('Data Loaded')

%% Load FOOOF Data - Group Probability

% Set filename
%filename = 'Group_Osc_Prob_Viz.mat';
filename = 'json_group_osc_prob_viz.mat';

% Load Group osc-prob data
load([datapath, filename])
fooof_dat = {[], theta_prob', alpha_prob', beta_prob'};
fooof_labels = {'xx', 'Theta_Prob', 'Alpha_Prob', 'Beta_Prob'};

disp('Data Loaded')

%% Load FOOOF Data - Osc Score

% Set filename
filename = 'json_group_osc_score_viz.mat';
%filename = 'diff_group_osc_score_viz.mat';

% Load Grop osc-score data
load([datapath, filename])
fooof_dat = {[], theta_score', alpha_score', beta_score'};
fooof_labels = {'xx', 'Theta_Score', 'Alpha_Score', 'Beta_Score'};

disp('Data Loaded')

%% Load FOOOF Data - Slopes

% Set filename
filename = 'Group_Slopes.mat';
%filename = 'Slopes_Age_Corr.mat';

% Load Slope data
load([datapath, filename]);

slopes = slopes';

fooof_dat = {slopes};
fooof_labels = {'Slopes'};

disp('Data Loaded');

%% Set FOOOF Data

for i = 1:length(fooof_dat)
    
    % Load BST dat file to use
    file_name = ['results_wMNE_MEG_KERNEL_160423_210', num2str(3+i), '.mat'];
    load([bst_path, file_name]);
    
    % Set data
    Time = 0;
    ImageGridAmp = fooof_dat{i};
    Comment = fooof_labels{i};

    % Save bst 
    save([bst_path, file_name], var_names{:});
    clear(var_names{:})

end

disp('Viz Set')

%%

% NOTE: SECTION BELOW ONLY HAD TO BE RUN ONCE TO CREATE A FILE TO USE
% DON'T RE-RUN UNLESS NEED THE BST VARS FILE AGAIN

%% Save a list of variable names to save a proper bst file
clear all; close all;
load('/Users/thomasdonoghue/Documents/brainstorm_db/Resting/data/Foof_vis/VIS/results_wMNE_MEG_160224_1414.mat');

% Get list of all variables, make cell array of bst variable names
vars = whos;
var_names = {};
for i = 1:length(vars)
    var_names{i} = vars(i).name;
end

% Save list of variable names
save('var_names.mat', 'var_names')

