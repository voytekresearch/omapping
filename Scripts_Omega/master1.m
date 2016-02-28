
sFiles = {...
    'Subject_981954/@rawMNI0001_MEGs0003_resting_20121128_01_AUX/data_0raw_MNI0001_MEGs0003_resting_20121128_01_AUX.mat'};
badchan=input('Enter bad channels like so- {''channel1'',''channel2''}   ');
bad_channels(sFiles,badchan);


sFiles = {...
   'Subject_369737/@rawMNI0007_MEGs0002_noise_20130315_06/data_0raw_MNI0007_MEGs0002_noise_20130315_06.mat'};
noise_continuous(sFiles);


sFiles = {...
    'Subject_369737/@rawMNI0007_MEGs0002_noise_20130315_06/data_0raw_MNI0007_MEGs0002_noise_20130315_06.mat', ...
    'Subject_369737/@rawMNI0007_MEGs0002_resting_20130315_01_AUX/data_0raw_MNI0007_MEGs0002_resting_20130315_01_AUX.mat'};
notch_filter(sFiles);


sFiles = {...
    'Subject_369737/@rawMNI0007_MEGs0002_noise_20130315_06/data_0raw_MNI0007_MEGs0002_noise_20130315_06.mat', ...
    'Subject_369737/@rawMNI0007_MEGs0002_resting_20130315_01_AUX/data_0raw_MNI0007_MEGs0002_resting_20130315_01_AUX.mat'};
delete_folder(sFiles);


sFiles = {...
    'Subject_369737/@rawMNI0007_MEGs0002_resting_20130315_01_AUX_notch/data_0raw_MNI0007_MEGs0002_resting_20130315_01_AUX_notch.mat'};
artifact_detection(sFiles);


