% Load PSD data and set to visualize

%% Set Up

% Load a list of variable names for a Bst file
load('var_names.mat')

% Set paths
datapath = '/Users/thomasdonoghue/Documents/Research/1-Projects/OMEGA/2-Data/MEG/4-Viz/';
bst_path = '/Users/thomasdonoghue/Documents/brainstorm_db/om_viz/data/om_viz_subj/Visualizing-MEG-FOOF/';

%% Load FOOF Data - Oscillation Frequencies

subj = 390845;
filename = [num2str(subj), '_Foof_Viz.mat'];

% Load FOOF Data
load([datapath, filename]);
%clear filename datapath

%
foof_dat = {slopes, thetas, alphas, betas, lowgammas};
foof_labels = {['Slopes_', num2str(subj)], ['Thetas_', num2str(subj)], ...
    ['Alphas_', num2str(subj)], ['Betas_', num2str(subj)], ['LowGammas_', num2str(subj)]};

disp('Data Loaded')

%% Load FOOF Data - Group Probability

% Set filename
filename = 'Group_Osc_Prob_Viz.mat';

% Load Group osc-prob data
load([datapath, filename])
foof_dat = {[], theta_prob, alpha_prob, beta_prob, lowgamma_prob};
foof_labels = {'xx', 'Theta_Prob', 'Alpha_Prob', 'Beta_Prob', 'LowGamma_Prob'};

disp('Data Loaded')

%% Set FOOF Data

for i = 1:length(foof_dat)
    
    % Load BST dat file to use
    file_name = ['results_wMNE_MEG_KERNEL_160423_210', num2str(3+i), '.mat'];
    load([bst_path, file_name]);
    
    % Set data
    Time = 0;
    ImageGridAmp = foof_dat{i};
    Comment = foof_labels{i};

    % Save bst 
    save([bst_path, file_name], var_names{:});
    clear(var_names{:})

end

disp('Viz Set')

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

