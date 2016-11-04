
SubjectNames = {...
    'Subject_981954'};
RawFilesanat = {...
    'C:\Users\pseb1015\Documents\Omega_Mapping\Subjects_Data\Subject_981954\Anatomy'};
RawFilesdata = {...
    'C:\Users\pseb1015\Documents\Omega_Mapping\Subjects_Data\Subject_369737\Data\OMEGA_V01_369737_MNI0007_MEGs0002_resting_20130315_01_AUX.ds_001\MNI0007_MEGs0002_resting_20130315_01_AUX.ds'};

nas = [123, 220, 112];
lpa = [43, 109, 107];
rpa = [212, 116, 98];
ac = [129, 133, 138];
pc = [131, 104, 137];
ih = [131, 114, 196];

% Import Anatomy

import_anatomy(SubjectNames,RawFilesanat,nas,lpa,rpa,ac,pc,ih);


% Link to raw file 

link_raw(SubjectNames,RawFilesdata);


% PSD - Channel Outlier

sFiles = {...
    'Subject_981954/@rawMNI0001_MEGs0003_resting_20121128_01_AUX/data_0raw_MNI0001_MEGs0003_resting_20121128_01_AUX.mat'};
PSD_COC();

copyfile('C:\Users\pseb1015\Documents\Omega_Mapping\brainstorm_db\Protocol01\data\Subject_981954\@rawMNI0001_MEGs0003_resting_20121128_01_AUX\timefreq_psd_160226_2028.mat', 'C:\Users\pseb1015\Documents\Omega_Mapping\Subjects_Processed\Subject_981954\psd_coc.mat')

