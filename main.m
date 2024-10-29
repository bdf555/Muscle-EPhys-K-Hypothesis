

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format shortG
%format long
addpath(pwd);
addpath(fullfile(pwd,'channel defs'));

clear('Chan_base');
clear('Molecs_base');
clear('protocol');
clear('current_stim_freq')
clear('add_capac_coupling')

paramfile = 'p_Rich_Kidea_19wt_norm_to_low_Struyk_v6b';
%paramfile = 'p_Rich_Kidea_19wt_norm_to_high_Struyk_v6b';
%paramfile = 'p_Rich_Kidea_19hypoK_norm_to_low_Struyk_v6b';
%paramfile = 'p_Rich_Kidea_19hypoK_norm_to_high_Struyk_v6b';
%paramfile = 'p_Rich_Kidea_19hyperK_norm_to_low_Struyk_v6b';
%paramfile = 'p_Rich_Kidea_19hyperK_norm_to_high_Struyk_v6b';
%paramfile = 'p_Rich_Kidea_19at_norm_to_low_Struyk_v6b';
%paramfile = 'p_Rich_Kidea_19at_norm_to_high_Struyk_v6b';


run(fullfile(pwd,'parameter files',paramfile));

process_params


if exist("filename",'var')
    load_full_data
else
    FileData = [];
end

fitData = format_fit_data_res_inc(FileData, uIpms);


    if do_fits_or_onerun == 1   % fit

        [all_params, outstruct] = fit_start(fitData);        

    elseif do_fits_or_onerun == 2   % one run

        if exist('onerun_parms','var')
            all_params = onerun_parms;
        else
            all_params = fitData.uIpms.fit.param_list;
        end
        outstruct = mep_runsim(all_params, fitData);
        outstruct.all_params = all_params;

    end


    if do_fits_or_onerun ~= 0
        if ~isfolder('model runs')
            mkdir('model runs');
        end
        svfilename = [save_file_prefix, '_', datestr(now,30)];
        save(fullfile(pwd,'model runs',svfilename),'all_params','uIpms','outstruct');
    end

    output_mep

rmpath(pwd);
rmpath(fullfile(pwd,'channel defs'));
