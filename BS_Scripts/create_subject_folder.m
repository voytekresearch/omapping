%% Create folder with subject data 

function create_subject_folder()
mkdir test_original;
for i= 1 : 10
    flname = strcat('test_original/Subject',int2str( i))
    mkdir (flname);
end
end

%% Move(or copy) subject into EMPTY subject folders
