% This script should be run after master2 and bad segments have been
% marked.

pipelinecount=4;

% Import recordings as 5s blocks
sFiles=inputfiles(pipelinecount);
import_meg(sFiles);
pipelinecount=pipelinecount+1;

% Remove DC offset
sFiles=inputfiles(pipelinecount);
dc_offset_block(sFiles);
pipelinecount=pipelinecount+1;

% Resample to 1000 Hz
sFiles=inputfiles(pipelinecount);
resample(sFiles);
pipelinecount=pipelinecount+1;

% Compute Noise Covariance
sFiles=inputfiles(pipelinecount);
noise_cov(sFiles);
pipelinecount=pipelinecount+1;

% Compute Head Model
sFiles=inputfiles(pipelinecount);
head_model(sFiles);
pipelinecount=pipelinecount+1;

% Compute Sources
sFiles=inputfiles(pipelinecount);
sources(sFiles);
pipelinecount=pipelinecount+1;

% Project on Default Anatomy
sFiles=inputfiles(pipelinecount);
project_da(sFiles);
pipelinecount=pipelinecount+1;

% Project on Default Anatomy
sFiles=inputfiles(pipelinecount);
psd_sensor(sFiles);
pipelinecount=pipelinecount+1;




sFiles = {...
    'Subject_981954/MNI0001_MEGs0003_resting_20121128_01_AUX_notch_bl/timefreq_psd_160227_0042.mat', ...
    'Subject_981954/MNI0001_MEGs0003_resting_20121128_01_AUX_notch_bl/timefreq_psd_160227_0043.mat', ...
    'Subject_981954/MNI0001_MEGs0003_resting_20121128_01_AUX_notch_bl/timefreq_psd_160227_0044.mat', ...
    'Subject_981954/MNI0001_MEGs0003_resting_20121128_01_AUX_notch_bl/timefreq_psd_160227_45.mat', ...
    'Subject_981954/MNI0001_MEGs0003_resting_20121128_01_AUX_notch_bl/timefreq_psd_160227_0045.mat', ...
    'Subject_981954/MNI0001_MEGs0003_resting_20121128_01_AUX_notch_bl/timefreq_psd_160227_0046.mat', ...
    'Subject_981954/MNI0001_MEGs0003_resting_20121128_01_AUX_notch_bl/timefreq_psd_160227_0047.mat', ...
    'Subject_981954/MNI0001_MEGs0003_resting_20121128_01_AUX_notch_bl/timefreq_psd_160227_48.mat', ...}
average_psd_files(sFiles);
delete_files(sFiles);


sFiles = {...
    'link|Subject_981954/MNI0001_MEGs0003_resting_20121128_01_AUX_notch_bl/results_wMNE_MEG_KERNEL_160227_0025.mat|Subject_981954/MNI0001_MEGs0003_resting_20121128_01_AUX_notch_bl/data_block002_bl_02.mat', ...
    'link|Subject_981954/MNI0001_MEGs0003_resting_20121128_01_AUX_notch_bl/results_wMNE_MEG_KERNEL_160227_0025.mat|Subject_981954/MNI0001_MEGs0003_resting_20121128_01_AUX_notch_bl/data_block005_bl_02.mat', ...
    'link|Subject_981954/MNI0001_MEGs0003_resting_20121128_01_AUX_notch_bl/results_wMNE_MEG_KERNEL_160227_0025.mat|Subject_981954/MNI0001_MEGs0003_resting_20121128_01_AUX_notch_bl/data_block006_bl_02.mat', ...
    'link|Subject_981954/MNI0001_MEGs0003_resting_20121128_01_AUX_notch_bl/results_wMNE_MEG_KERNEL_160227_0025.mat|Subject_981954/MNI0001_MEGs0003_resting_20121128_01_AUX_notch_bl/data_block007_bl_02.mat', ...
    'link|Subject_981954/MNI0001_MEGs0003_resting_20121128_01_AUX_notch_bl/results_wMNE_MEG_KERNEL_160227_0025.mat|Subject_981954/MNI0001_MEGs0003_resting_20121128_01_AUX_notch_bl/data_block011_bl_02.mat', ...
    'link|Subject_981954/MNI0001_MEGs0003_resting_20121128_01_AUX_notch_bl/results_wMNE_MEG_KERNEL_160227_0025.mat|Subject_981954/MNI0001_MEGs0003_resting_20121128_01_AUX_notch_bl/data_block012_bl_02.mat', ...
    'link|Subject_981954/MNI0001_MEGs0003_resting_20121128_01_AUX_notch_bl/results_wMNE_MEG_KERNEL_160227_0025.mat|Subject_981954/MNI0001_MEGs0003_resting_20121128_01_AUX_notch_bl/data_block013_bl_02.mat', ...
    'link|Subject_981954/MNI0001_MEGs0003_resting_20121128_01_AUX_notch_bl/results_wMNE_MEG_KERNEL_160227_0025.mat|Subject_981954/MNI0001_MEGs0003_resting_20121128_01_AUX_notch_bl/data_block014_bl_02.mat', ...
    'link|Subject_981954/MNI0001_MEGs0003_resting_20121128_01_AUX_notch_bl/results_wMNE_MEG_KERNEL_160227_0025.mat|Subject_981954/MNI0001_MEGs0003_resting_20121128_01_AUX_notch_bl/data_block017_bl_02.mat', ...
    'link|Subject_981954/MNI0001_MEGs0003_resting_20121128_01_AUX_notch_bl/results_wMNE_MEG_KERNEL_160227_0025.mat|Subject_981954/MNI0001_MEGs0003_resting_20121128_01_AUX_notch_bl/data_block018_bl_02.mat', ...
    'link|Subject_981954/MNI0001_MEGs0003_resting_20121128_01_AUX_notch_bl/results_wMNE_MEG_KERNEL_160227_0025.mat|Subject_981954/MNI0001_MEGs0003_resting_20121128_01_AUX_notch_bl/data_block019_bl_02.mat', ...}
PSD_source(sFiles);












