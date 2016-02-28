% After bad segments marked

sFiles = {...
    'Subject_981954/@rawMNI0001_MEGs0003_resting_20121128_01_AUX_notch/data_0raw_MNI0001_MEGs0003_resting_20121128_01_AUX_notch.mat'};
dc_offset(sFiles);


sFiles = {...
    'Subject_981954/@rawMNI0001_MEGs0003_resting_20121128_01_AUX_notch_bl/data_0raw_MNI0001_MEGs0003_resting_20121128_01_AUX_notch_bl.mat'};
sensor_psd(sFiles);


sFiles = {...
    'Subject_981954/@rawMNI0001_MEGs0003_resting_20121128_01_AUX_notch_bl/data_0raw_MNI0001_MEGs0003_resting_20121128_01_AUX_notch_bl.mat'};
SubjectNames = {...
    'Subject_981954'};
import_meg(sFiles,SubjectNames);


sFiles = {...
    'Subject_981954/MNI0001_MEGs0003_resting_20121128_01_AUX_notch_bl/data_block002.mat', ...
    'Subject_981954/MNI0001_MEGs0003_resting_20121128_01_AUX_notch_bl/data_block005.mat', ...
    'Subject_981954/MNI0001_MEGs0003_resting_20121128_01_AUX_notch_bl/data_block006.mat', ...
    'Subject_981954/MNI0001_MEGs0003_resting_20121128_01_AUX_notch_bl/data_block007.mat', ...
    'Subject_981954/MNI0001_MEGs0003_resting_20121128_01_AUX_notch_bl/data_block011.mat', ...
    'Subject_981954/MNI0001_MEGs0003_resting_20121128_01_AUX_notch_bl/data_block012.mat', ...
    'Subject_981954/MNI0001_MEGs0003_resting_20121128_01_AUX_notch_bl/data_block058.mat', ...
    'Subject_981954/MNI0001_MEGs0003_resting_20121128_01_AUX_notch_bl/data_block059.mat'};

dc_offset_block(sFiles);





sFiles = {...
    'Subject_981954/@rawMNI0001_MEGs0003_noise_20121128_02_notch/data_0raw_MNI0001_MEGs0003_noise_20121128_02_notch.mat'};
noise_cov(sFiles);


sFiles = {...
    'Subject_981954/MNI0001_MEGs0003_resting_20121128_01_AUX_notch_bl/data_block002_bl_02.mat', ...
    'Subject_981954/MNI0001_MEGs0003_resting_20121128_01_AUX_notch_bl/data_block005_bl_02.mat', ...
    'Subject_981954/MNI0001_MEGs0003_resting_20121128_01_AUX_notch_bl/data_block006_bl_02.mat', ...
    'Subject_981954/MNI0001_MEGs0003_resting_20121128_01_AUX_notch_bl/data_block007_bl_02.mat', ...
    'Subject_981954/MNI0001_MEGs0003_resting_20121128_01_AUX_notch_bl/data_block011_bl_02.mat', ...
    'Subject_981954/MNI0001_MEGs0003_resting_20121128_01_AUX_notch_bl/data_block012_bl_02.mat', ...
    'Subject_981954/MNI0001_MEGs0003_resting_20121128_01_AUX_notch_bl/data_block013_bl_02.mat', ...
    'Subject_981954/MNI0001_MEGs0003_resting_20121128_01_AUX_notch_bl/data_block039_bl_02.mat', ...
    'Subject_981954/MNI0001_MEGs0003_resting_20121128_01_AUX_notch_bl/data_block040_bl_02.mat', ...
    'Subject_981954/MNI0001_MEGs0003_resting_20121128_01_AUX_notch_bl/data_block041_bl_02.mat', ...
    'Subject_981954/MNI0001_MEGs0003_resting_20121128_01_AUX_notch_bl/data_block042_bl_02.mat', ...
    'Subject_981954/MNI0001_MEGs0003_resting_20121128_01_AUX_notch_bl/data_block059_bl_02.mat'};
head(sFiles);


sFiles = {...
    'link|Subject_369737/MNI0007_MEGs0002_resting_20130315_01_AUX_notch_bl/results_wMNE_MEG_KERNEL_160227_1643.mat|Subject_369737/MNI0007_MEGs0002_resting_20130315_01_AUX_notch_bl/data_block002_bl.mat', ...
    'link|Subject_369737/MNI0007_MEGs0002_resting_20130315_01_AUX_notch_bl/results_wMNE_MEG_KERNEL_160227_1643.mat|Subject_369737/MNI0007_MEGs0002_resting_20130315_01_AUX_notch_bl/data_block003_bl.mat', ...
    'link|Subject_369737/MNI0007_MEGs0002_resting_20130315_01_AUX_notch_bl/results_wMNE_MEG_KERNEL_160227_1643.mat|Subject_369737/MNI0007_MEGs0002_resting_20130315_01_AUX_notch_bl/data_block004_bl.mat', ...
    'link|Subject_369737/MNI0007_MEGs0002_resting_20130315_01_AUX_notch_bl/results_wMNE_MEG_KERNEL_160227_1643.mat|Subject_369737/MNI0007_MEGs0002_resting_20130315_01_AUX_notch_bl/data_block005_bl.mat', ...
    'link|Subject_369737/MNI0007_MEGs0002_resting_20130315_01_AUX_notch_bl/results_wMNE_MEG_KERNEL_160227_1643.mat|Subject_369737/MNI0007_MEGs0002_resting_20130315_01_AUX_notch_bl/data_block008_bl.mat', ...
    'link|Subject_369737/MNI0007_MEGs0002_resting_20130315_01_AUX_notch_bl/results_wMNE_MEG_KERNEL_160227_1643.mat|Subject_369737/MNI0007_MEGs0002_resting_20130315_01_AUX_notch_bl/data_block009_bl.mat', ...
    'link|Subject_369737/MNI0007_MEGs0002_resting_20130315_01_AUX_notch_bl/results_wMNE_MEG_KERNEL_160227_1643.mat|Subject_369737/MNI0007_MEGs0002_resting_20130315_01_AUX_notch_bl/data_block011_bl.mat', ...
    'link|Subject_369737/MNI0007_MEGs0002_resting_20130315_01_AUX_notch_bl/results_wMNE_MEG_KERNEL_160227_1643.mat|Subject_369737/MNI0007_MEGs0002_resting_20130315_01_AUX_notch_bl/data_block012_bl.mat', ...
    'link|Subject_369737/MNI0007_MEGs0002_resting_20130315_01_AUX_notch_bl/results_wMNE_MEG_KERNEL_160227_1643.mat|Subject_369737/MNI0007_MEGs0002_resting_20130315_01_AUX_notch_bl/data_block013_bl.mat', ...}
power_per_timeblock(sFiles);


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