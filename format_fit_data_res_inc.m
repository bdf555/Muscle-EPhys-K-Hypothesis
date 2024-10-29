function fitData = format_fit_data_res_inc(fullData, uIpms)

    fitData.uIpms = uIpms;
    
    % list of time-dependent variables
    % number of elements for each, and the order they appear in vectors
    
    % elements are:
    %       1st: total number of instances
    %       2nd: first index in long vector
    %       3rd: last index in long vector
    %       4th: number of segments (1 or ns)
    %       5th: number of shell instances (1 or nsh)
    
    num_segments = uIpms.gen.num_segments;
    num_shells = uIpms.fit.num_shells;
    
        li = 0;
    fitData.tdvars{1} = [num_segments, li+1, li+num_segments, num_segments, 1];     % Vm
        li = li + fitData.tdvars{1}(1);
    fitData.tdvars{2} = [1, li+1, li+1, 1, 1];     % V0vps
        li = li + fitData.tdvars{2}(1);
    fitData.tdvars{3} = [1, li+1, li+1, 1, 1];     % dV0vps
        li = li + fitData.tdvars{3}(1);
    fitData.tdvars{4} = [num_segments*num_shells, li+1, li+num_segments*num_shells, num_segments, num_shells];     % Vt
        li = li + fitData.tdvars{4}(1);
    fitData.tdvars{5} = [num_segments*num_shells, li+1, li+num_segments*num_shells, num_segments, num_shells];     % Volume of t-tubules
        li = li + fitData.tdvars{5}(1);
    fitData.tdvars{6} = [num_segments, li+1, li+num_segments, num_segments, 1];     % Volume of intracellular compartments
        li = li + fitData.tdvars{6}(1);
    fitData.tdvars{7} = [1, li+1, li+1, 1, 1];     % Volume of extracellular space
        li = li + fitData.tdvars{7}(1);
    state_num = 7;   
        
    fitData.num_molecs = length(uIpms.Molecs_base);
    for nmolecvs = 1:fitData.num_molecs
        state_num = state_num + 1;       
        fitData.tdvars{state_num} = [num_segments*num_shells, li+1, li+num_segments*num_shells, num_segments, num_shells];     % t conc
        li = li + fitData.tdvars{state_num}(1);
        
        state_num = state_num + 1;       
        fitData.tdvars{state_num} = [num_segments, li+1, li+num_segments, num_segments, 1];     % i conc
        li = li + fitData.tdvars{state_num}(1);
        
        state_num = state_num + 1;       
        fitData.tdvars{state_num} = [1, li+1, li+1, 1, 1];     % e conc
        li = li + fitData.tdvars{state_num}(1);
    end
    
        
    for nchs = 1:uIpms.num_chans
        num_states = uIpms.Chan_base(nchs).vars{1};
        for ins = 1:num_states
            state_num = state_num + 1;
            fitData.tdvars{state_num} = [num_segments, li+1, li+num_segments, num_segments, 1];     % sarc membrane state variable for channel
            li = li + fitData.tdvars{state_num}(1);
            
            state_num = state_num + 1;    
            fitData.tdvars{state_num} = [num_segments*num_shells, li+1, li+num_segments*num_shells, num_segments, num_shells];     % tt membrane state variable for channel
            li = li + fitData.tdvars{state_num}(1);
        end
    end
    
    fitData.tdvars_data = [length(fitData.tdvars) li];  % [number of distinct variables, total number of state variables (so multiplied by segments and shells)]
    
    R = 8.3144598;			% J/(mol*K)
    F = 96485.33;			% coulombs/mol
    fitData.RTF = (R*uIpms.gen.Temp/F)*1000;          % multiplied by 1000 to convert to mV
    RTF = fitData.RTF;
    TFactor = 1./ (3.^((37-(uIpms.gen.Temp - 273))./10));
    fitData.F = F;			% coulombs/mol
    fitData.H2O_molar_volume = 18;  % cm^3/mole
    
    fitData.deltax = uIpms.gen.cell_length ./ fitData.uIpms.gen.num_segments;
    fitData.Cm_to_C = 2.*pi.*uIpms.gen.cell_radius.*fitData.deltax .* 1000;  % converts to nF for membrane of one segment
    fitData.G_to_g = fitData.deltax .* 2 .* pi .* uIpms.gen.cell_radius .*1000;   % converts to uS of 1 segment cell membrane
    fitData.sarc_seg_area = fitData.deltax .* 2 .* pi .* uIpms.gen.cell_radius; % cm^2 of sarcolemma surface area
    fitData.R_to_r = 1./(fitData.deltax .* 2 .* pi .* uIpms.gen.cell_radius .*1000);   % converts kOhms cm^2 of Ra to MOhms

    if uIpms.fit.num_shells > 0
        fitData.delta_r = uIpms.gen.cell_radius/uIpms.fit.num_shells;
        radius_prev = 0;
        fitData.tt_total_surface_area = 0;
        fitData.radius_tt_compartment = zeros(uIpms.fit.num_shells, 1);
        fitData.volume_tt_compartment = zeros(uIpms.fit.num_shells, 1);
        fitData.area_longitudinal_anulus = zeros(uIpms.fit.num_shells, 1);
        fitData.surface_area_tt_compartment = zeros(uIpms.fit.num_shells, 1);	% membrane surface area between one t membrane subcompartment and the intracellular compartment for one fiber cable segment
        fitData.diff_cross_section_area_tt_compartment = zeros(uIpms.fit.num_shells,1);

        for j = 1:fitData.uIpms.fit.num_shells
            fitData.radius_tt_compartment(j) = fitData.delta_r * j;
            fitData.volume_tt_compartment(j) = uIpms.fit.rho_t * pi * fitData.deltax *(fitData.radius_tt_compartment(j)^2 - radius_prev^2);
            fitData.area_longitudinal_anulus(j) = pi * (fitData.radius_tt_compartment(j)^2 - radius_prev^2);
            fitData.surface_area_tt_compartment(j) = (fitData.volume_tt_compartment(j)) / uIpms.fit.zeta;
            fitData.tt_total_surface_area = fitData.tt_total_surface_area + fitData.surface_area_tt_compartment(j);
            radius_prev = radius_prev + fitData.delta_r; 
        end
        fitData.tt_total_volume = sum(uIpms.gen.num_segments.*fitData.volume_tt_compartment);    % total, over all segments
        fitData.intra_volume = pi .* uIpms.gen.cell_radius.^2 .* uIpms.gen.cell_length - sum(uIpms.gen.num_segments.*fitData.volume_tt_compartment);
        fitData.sarc_surf_area = 2.*pi.*uIpms.gen.cell_radius.*uIpms.gen.cell_length;
        fitData.tot_surf_area = fitData.sarc_surf_area + uIpms.gen.num_segments.*fitData.tt_total_surface_area;
    else
        fitData.intra_volume = pi .* uIpms.gen.cell_radius.^2 .* uIpms.gen.cell_length;
        fitData.sarc_surf_area = 2.*pi.*uIpms.gen.cell_radius.*uIpms.gen.cell_length;
        fitData.tot_surf_area = fitData.sarc_surf_area;
    end
    
    fitData.seg_intra_volume = fitData.intra_volume ./ uIpms.gen.num_segments;
    fitData.extra_volume = fitData.intra_volume .* uIpms.gen.extrac_volume_multiple;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    


    fitData.num_protocols = length(uIpms.protocols);
    fitData.prot = uIpms.protocols;
    
    fitData.traces_to_fit = [];

    for nprot = 1:fitData.num_protocols
        
        filenum = fitData.prot{nprot}.filenum;
        if filenum ~= 0
            thisfile = fullData{filenum};
        end
        
        fitData.prot{nprot}.chan_perm_factors = set_chan_perm_factors(fitData.prot{nprot}, fitData);
        
        if fitData.prot{nprot}.filenum == 0 && ismember(fitData.prot{nprot}.input,[1 4])
            error('can not use input methods that require data when doing pure sim')
        end

        if fitData.prot{nprot}.input == 1
            fitData.prot{nprot}.forced_time_evals = [];
        end

        if fitData.prot{nprot}.input == 2

            if 1
            ic_times = fitData.prot{nprot}.input_current_density(:,1);
            ic_init_val = fitData.prot{nprot}.input_current_density_start*2*pi*uIpms.gen.cell_radius*uIpms.gen.cell_length;
            ic_vals = zeros(length(ic_times),1);
            for nis = 1:length(ic_times)
                ic_vals(nis) = fitData.prot{nprot}.input_current_density(nis,2)*2*pi*uIpms.gen.cell_radius*uIpms.gen.cell_length;
            end
            
            ic_v = ic_init_val;
            old_ic_val = ic_init_val;
            forced_time_evals = [];
            
            for nic = 1:length(ic_times)
                idx_val = (nic-1)*2+2;
                ic_t(idx_val) = ic_times(nic) - 1e-10;
                ic_v(idx_val) = old_ic_val;
                ic_t(idx_val+1) = ic_times(nic);
                ic_v(idx_val+1) = ic_vals(nic);
                old_ic_val = ic_vals(nic);
                forced_time_evals = [forced_time_evals ic_times(nic) + 1e-7];
            end
            
            if fitData.prot{nprot}.filenum ~= 0
                t_sim_end = find( thisfile.times >= fitData.prot{nprot}.sim_times(2), 1 );		% index
                final_t = thisfile.times(t_sim_end);
            else
                final_t = fitData.prot{nprot}.final_t;
            end

            ic_t = [ic_t final_t];
            ic_v = [ic_v old_ic_val];
            
            fitData.prot{nprot}.ic_mat = [ic_t' ic_v'];
            
            % if events are short duration, need to force ode solver to evaluate during event.
            % typically, include time points in middle of current pulse
            %forced_time_evals = [50.1 75 100.1];
            fitData.prot{nprot}.forced_time_evals = forced_time_evals;
            end

        end
        
        if isfield(fitData.prot{nprot},'molecs_change')
            % construct interpolation table for changes to extracellular
            % ion concentrations
            for im = fitData.prot{nprot}.molecs_change.modified_molecs
                mc_t = [];
                mc_dmdt = [];
                mmat = fitData.prot{nprot}.molecs_change.molec(im).mat;
                mctimes = mmat(:,1);
                mcvals = mmat(:,2);
                mcdur = mmat(:,3);
                mc_t(1) = 0;
                mc_dmdt(1) = 0;
                for nim = 1:size(mmat,1)
                    dmdt = mcvals(nim).*fitData.extra_volume.*1e3./mcdur(nim);     % nanomoles/ms when mcvals is in mM, volume is in cm^3, mcdur is in ms

                    idx_val = (nim-1)*4+2;
                    mc_t(idx_val) = mctimes(nim) - 1e-7;
                    mc_dmdt(idx_val) = 0;
                    mc_t(idx_val+1) = mctimes(nim);
                    mc_dmdt(idx_val+1) = dmdt;
                    fitData.prot{nprot}.forced_time_evals = union(fitData.prot{nprot}.forced_time_evals, mctimes(nim) + 1e-7);

                    mc_t(idx_val+2) = mctimes(nim) + mcdur(nim) - 1e-7;
                    mc_dmdt(idx_val+2) = dmdt;
                    mc_t(idx_val+3) = mctimes(nim) + mcdur(nim);
                    mc_dmdt(idx_val+3) = 0;
                    fitData.prot{nprot}.forced_time_evals = union(fitData.prot{nprot}.forced_time_evals, mctimes(nim) + mcdur(nim) + 1e-7);
                    
                end

                if fitData.prot{nprot}.filenum ~= 0
                    t_sim_end = find( thisfile.times >= fitData.prot{nprot}.sim_times(2), 1 );		% index
                    final_t = thisfile.times(t_sim_end);
                else
                    final_t = fitData.prot{nprot}.final_t;
                end

                mc_t = [mc_t final_t];
                mc_dmdt = [mc_dmdt 0];

                fitData.prot{nprot}.mc_mat{im} = [mc_t' mc_dmdt'];
            end
        end



        if isfield(fitData.prot{nprot},'sim_voltage_input')
            simv_final_row = [fitData.prot{nprot}.sim_times(2) fitData.prot{nprot}.sim_voltage_input(end,2:end)];
            fitData.prot{nprot}.sim_voltage_input = [fitData.prot{nprot}.sim_voltage_input; simv_final_row];
        end
        
        if fitData.prot{nprot}.input == 6
            fitData.prot{nprot}.times_for_interp = fitData.prot{nprot}.sim_voltage_input(:,1);
            fitData.prot{nprot}.volt_traces_for_interp = fitData.prot{nprot}.sim_voltage_input(:,2:end);
        end
        
        if ~isfield(fitData.prot{nprot},'tshift')
            fitData.prot{nprot}.tshift = 0.0;
        end
        
        if ~isfield(fitData.prot{nprot},'steps_to_run')
            fitData.prot{nprot}.steps_to_run = 1;
            fitData.prot{nprot}.uniform_initial_steady_state = 1;
        end
        
        if ~isfield(fitData.prot{nprot},'time_to_establish_steady_state')
            fitData.prot{nprot}.time_to_establish_steady_state = 1e4;
        end
        
        if ~isfield(fitData.prot{nprot},'Vminit_method')
            if fitData.prot{nprot}.filenum ~= 0
                fitData.prot{nprot}.Vminit_method = 2;
            else
                fitData.prot{nprot}.Vminit_method = 3;
            end
        elseif fitData.prot{nprot}.Vminit_method == 1 && ~isfield(fitData.prot{nprot},'RM_orig')
            error('must specify RM_orig value when using Vminit_method = 1')
        elseif fitData.prot{nprot}.Vminit_method == 2 && fitData.prot{nprot}.filenum == 0
            error('Impossible Vminit_method.  Cannot use zeros from experimental trace as Vminit when doing pure sim')
        elseif fitData.prot{nprot}.Vminit_method == 4 && ~ismember(fitData.prot{nprot}.input,[3,4,5,6])
            error('Invalid Vminit_method.  Cannot base Vminit on generated V if input not from generated V')
        end

        if fitData.prot{nprot}.filenum ~= 0
        
            t_sim_start = find( thisfile.times >= fitData.prot{nprot}.sim_times(1), 1);		% index
            t_sim_end = find( thisfile.times >= fitData.prot{nprot}.sim_times(2), 1 );		% index
            t_sim_idxs = [t_sim_start: t_sim_end];
            t_interp_idxs = [t_sim_start - 10 : t_sim_end + 10];

            res_inc = uIpms.fit.time_resolution_increase;   % for adding more time points
            delta_t = thisfile.times(2) - thisfile.times(1);
            new_delta_t = delta_t / res_inc;
            times_to_sim_temp = thisfile.times(t_sim_idxs);
            
            fitData.prot{nprot}.times_to_sim = [times_to_sim_temp(1) : new_delta_t : times_to_sim_temp(end)]';

            t_fit_idxs = [];
            for j = 1:size(fitData.prot{nprot}.fit_times,1)
                t_fit_start = find( thisfile.times >= fitData.prot{nprot}.fit_times(j,1), 1);		% index of tfit_start
                t_fit_end = find( thisfile.times >= fitData.prot{nprot}.fit_times(j,2), 1 );		% index of tfit_end
                t_fit_idxs = [t_fit_idxs t_fit_start : t_fit_end];
            end
            times_to_fit_temp = thisfile.times(t_fit_idxs);
            [~, fitData.prot{nprot}.t_fit_idxs_of_sim_times] = ismember(t_fit_idxs,t_sim_idxs);

            fitData.prot{nprot}.times_to_fit = thisfile.times(t_fit_idxs);

            zero_region_start = find( thisfile.times >= fitData.prot{nprot}.zero_times_before(1), 1);		% index of zero region start
            zero_region_end = find( thisfile.times >= fitData.prot{nprot}.zero_times_before(2), 1 );		% index of zero region end
            zero_idxs = [zero_region_start: zero_region_end];

            
            volt_traces_to_sim_temp = thisfile.voltage_traces(t_sim_idxs, fitData.prot{nprot}.steps_to_run);
            fitData.prot{nprot}.volt_traces_to_sim = interp1(times_to_sim_temp, volt_traces_to_sim_temp, fitData.prot{nprot}.times_to_sim);

            if fitData.prot{nprot}.input ~= 6
                fitData.prot{nprot}.times_for_interp = thisfile.times(t_interp_idxs);
                fitData.prot{nprot}.volt_traces_for_interp = thisfile.voltage_traces(t_interp_idxs, fitData.prot{nprot}.steps_to_run);
            end
            
            volt_traces_to_fit_temp = thisfile.voltage_traces(t_fit_idxs, fitData.prot{nprot}.steps_to_run);
            fitData.prot{nprot}.volt_traces_to_fit = interp1(times_to_fit_temp, volt_traces_to_fit_temp, fitData.prot{nprot}.times_to_fit);
            
            if isfield(fitData.prot{nprot},'zero_times_after')
                zero_region_start = find( thisfile.times >= fitData.prot{nprot}.zero_times_after(1), 1);		% index of zero region start
                zero_region_end = find( thisfile.times >= fitData.prot{nprot}.zero_times_after(2), 1 );		% index of zero region end
                zero_final_idxs = [zero_region_start: zero_region_end];

                mean_final_current = mean(thisfile.current_traces(zero_final_idxs, fitData.prot{nprot}.steps_to_run));
            else
                mean_final_current = 0;
            end
            
            fitData.prot{nprot}.currents_to_fit = thisfile.current_traces(t_fit_idxs, fitData.prot{nprot}.steps_to_run) - mean_final_current;
            currents_to_sim_temp = thisfile.current_traces(t_sim_idxs, fitData.prot{nprot}.steps_to_run) - mean_final_current;       % only used for plotting
            fitData.prot{nprot}.currents_to_sim = interp1(times_to_sim_temp, currents_to_sim_temp, fitData.prot{nprot}.times_to_sim);       % only used for plotting
            fitData.prot{nprot}.currents_for_interp = thisfile.current_traces(t_interp_idxs, fitData.prot{nprot}.steps_to_run) - mean_final_current;
            fitData.prot{nprot}.mean_initial_current = mean(thisfile.current_traces(zero_idxs, fitData.prot{nprot}.steps_to_run));

            if ismember(fitData.prot{nprot}.input, [0 3 4 5 6])
                % voltage is input, current is output
                fitData.prot{nprot}.traces_to_fit = fitData.prot{nprot}.currents_to_fit(:);
                
                if isfield(fitData.prot{nprot}, 'step_label_time')
                    label_idx = find( thisfile.times >= fitData.prot{nprot}.step_label_time, 1);
                    fitData.prot{nprot}.voltage_steps = thisfile.voltage_traces(label_idx,:);
                end
                % fitData.prot{nprot}.steps_to_fit is really step_voltage values
                % it is only used for sigmoid steps via sP.Vc_mag in input_method = 4
                
            elseif fitData.prot{nprot}.input == 1 || fitData.prot{nprot}.input == 2
                % input is current, output is voltage
                fitData.prot{nprot}.traces_to_fit = fitData.prot{nprot}.volt_traces_to_fit(:);
            end
            
            
        else    % fitData.prot{nprot}.filenum = 0, pure sim no experimental data
            
            fitData.prot{nprot}.times_to_sim = [0 : fitData.prot{nprot}.created_delta_t : fitData.prot{nprot}.final_t]';
            if fitData.prot{nprot}.input == 6
                fitData.prot{nprot}.mean_initial_current = 0;   % might need to create new parameter for this, or calculate it based on current needed to "force" chosen initial voltage.
            elseif fitData.prot{nprot}.input == 2
                fitData.prot{nprot}.mean_initial_current = fitData.prot{nprot}.ic_mat(1,2);
            elseif fitData.prot{nprot}.input == 5 || fitData.prot{nprot}.input == 3
                fitData.prot{nprot}.mean_initial_current = 0;   % might need to create new parameter for this, or calculate it based on current needed to "force" chosen initial voltage.
            elseif fitData.prot{nprot}.input == 7
                fitData.prot{nprot}.mean_initial_current = fitData.prot{nprot}.ic_mat(1,2);
                
            end
            
        end
        
        if fitData.prot{nprot}.Vminit_method == 1
            fitData.prot{nprot}.mean_initial_voltage = fitData.prot{nprot}.RM_orig;
        elseif fitData.prot{nprot}.Vminit_method == 2
            fitData.prot{nprot}.mean_initial_voltage = mean(thisfile.voltage_traces(zero_idxs, fitData.prot{nprot}.steps_to_run));
        elseif fitData.prot{nprot}.Vminit_method == 3
            fitData.prot{nprot}.mean_initial_voltage = fitData.RTF*log(uIpms.Molecs_base(2).init_extrac ./ uIpms.Molecs_base(2).init_intrac) ./ uIpms.Molecs_base(2).valence;
        elseif fitData.prot{nprot}.Vminit_method == 4
            if fitData.prot{nprot}.input == 5 || fitData.prot{nprot}.input == 3
                fitData.prot{nprot}.mean_initial_voltage = fitData.prot{nprot}.sig_step_init;
            elseif fitData.prot{nprot}.input == 6 
                fitData.prot{nprot}.mean_initial_voltage = fitData.prot{nprot}.sim_voltage_input(1,2:end);
            end
        end
        
        fitData.prot{nprot} = init_concs(fitData.prot{nprot}, fitData);
        
        
    end
    

    
    
