
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_folder = fullfile(pwd,'experimental data');

FileData = [];

for nf = 1:length(uIpms.filenames)
    
    if isfile(fullfile(data_folder, [uIpms.filenames{nf} '.mat']))
        load(fullfile(data_folder, uIpms.filenames{nf}));
        FileData{nf} = fullData;
    else
        error('file does not exist');
    end
end

