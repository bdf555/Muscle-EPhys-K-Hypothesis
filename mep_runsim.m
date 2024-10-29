
function outstruct = mep_runsim(all_params, fitData)

    % sP stands for sim paramters.  All the numbers which don't change
    % during a single time-course simulation of the membrane.

    ns = fitData.uIpms.gen.num_segments;
    nsh = fitData.uIpms.fit.num_shells;
    
    sP.try_small_matrix = 1;
    
    sP.use_IC_comps = fitData.uIpms.fit.use_IC_comps;
    sP.track_volumes = fitData.uIpms.track_volumes;
    sP.track_icec_ion_numbers = fitData.uIpms.track_icec_ion_numbers;
    sP.track_tt_ion_numbers = fitData.uIpms.track_tt_ion_numbers;
    sP.use_electrical_drift = fitData.uIpms.use_electrical_drift;
    
    sP.RTF = fitData.RTF;
    sP.F = fitData.F;
    sP.num_molecs = fitData.num_molecs;
    sP.molec_chan = fitData.uIpms.molec_chan;
    
    sP.gen_params = all_params(fitData.uIpms.gen_param_idxs);
    
    sP.chan_params = [];
    if fitData.uIpms.each_step_diff_exp == 0
        for ncs = 1:fitData.uIpms.num_chans
            sP.chan_params{ncs} = all_params(fitData.uIpms.chan_param_idxs{ncs});
        end
    else
        for ncs = 1:fitData.uIpms.num_chans
            sP.chan_params{ncs}(1) = all_params(fitData.uIpms.fit.param_list_idxnum{ncs}{1}{1}{1}); % only 1 Na permeability value
        end
        for apidx = 6:length(all_params)
            prmdta = fitData.uIpms.fit.param_locator(apidx,:);   %[nCs nprots npns nstepsnum]
            if prmdta(3) ~= 1   % if not the permeability parameter
                sP.chan_params_vc{prmdta(2)}{prmdta(1)}(prmdta(3),prmdta(4)) = all_params(apidx);   % each row is a parameter, each column is a step
            end
        end
    end
    sP.chan_data = fitData.uIpms.Chan_base;
    sP.num_chans = fitData.uIpms.num_chans;


    sP.C = sP.gen_params(1) .* fitData.Cm_to_C;    % C per segment, nF
    sP.g_S = sP.gen_params(2) .* fitData.sarc_seg_area .* 1000;   % series conductance per segment, uS
    sP.H2O_perm_factor = sP.gen_params(5) .* fitData.sarc_seg_area .* fitData.H2O_molar_volume / 10^9;    % (Liters * cm^3) / (millimoles * ms)
                                                                                                       % multiplying by mM (osmolarity) gives cm^3/ms of water flux
    if nsh > 0 && sP.use_IC_comps
        sP.gIL = (fitData.area_longitudinal_anulus .* sP.gen_params(3) .* 1000) ./ fitData.deltax; % conductance (uS) between adjacent intracellular longitudinal segments
    else
        sP.gIL = sP.gen_params(3) .* 1000 .* pi .* fitData.uIpms.gen.cell_radius.^2 ./ fitData.deltax;
    end

    sP.chan_partial_cond_factor = [];
    for ncs = 1:fitData.uIpms.num_chans
        sP.chan_partial_cond_factor(ncs) = sP.chan_params{ncs}(1) .* fitData.sarc_seg_area .* 1000; % becomes uS/mM per segment if use_GHK, else uS per segment in find_max_cond_factor.m
                                                                                                    % becomes nA for channels where the max conductance (in find_max_cond_factor) doesn't
                                                                                                    % depend on molecule concentrations.
    end
        
    if nsh > 0
        
        G_bar_L = fitData.uIpms.fit.rho_t * fitData.uIpms.fit.sigma_t * sP.gen_params(4);	% mS/cm, t-tubular cable conductivity
        sP.gT = (2 .* pi .* fitData.radius_tt_compartment .* fitData.deltax .* G_bar_L .* 1000) ./ fitData.delta_r; % conductance (uS) between adjacent t tubule radial segments
        sP.gTL = (fitData.area_longitudinal_anulus .* G_bar_L .* 1000) ./ fitData.deltax; % conductance (uS) between adjacent t tubule longitudinal segments
        
        sP.gIR = (2 .* pi .* fitData.radius_tt_compartment .* fitData.deltax .* sP.gen_params(3) .* 1000) ./ fitData.delta_r; % conductance (uS) between adjacent intracellular radial segments
        
        sP.chant_partial_cond_factor = [];
        for ncs = 1:fitData.uIpms.num_chans
            sP.chant_partial_cond_factor{ncs} = sP.chan_params{ncs}(1) .* sP.chan_data(ncs).eta .* ones(ns,nsh) .* fitData.surface_area_tt_compartment' .* 1000;    % max conductance of each t-tubule compartment in 1 longitudinal segment. (uS/mM) if use_GHK, else uS
                                                                                                    % becomes nA for channels where the max conductance (in find_max_cond_factor) doesn't
                                                                                                    % depend on molecule concentrations.
            
        end
        
        sP.C_t = sP.gen_params(1) .* ones(ns,nsh) .* fitData.surface_area_tt_compartment' .* 1000;  % membrane capacitiance of each t-tubule compartment in 1 segment (nF)
        sP.C_t_seg = sP.gen_params(1) .* sum(fitData.surface_area_tt_compartment) .* 1000;      % total membrane capacitance of t-tubule compartment (all shells) for 1 segment (nF)
        if fitData.uIpms.fit.use_access_resistance
            sP.ra = fitData.uIpms.fit.Ra .* fitData.R_to_r;   % access resistance for 1 segment (MOhms)
        else
            sP.ra = 1./(2.*sP.gT(end));   % access resistance for 1 segment (MOhms)
        end
        
        sP.radius_tt_compartment = fitData.radius_tt_compartment;
        sP.radius_tt_for_diff = [0;sP.radius_tt_compartment];
        sP.radius_tt_for_diff = repmat(sP.radius_tt_for_diff,1,fitData.uIpms.gen.num_segments)';
        
        sP.volume_tt_compartment = fitData.volume_tt_compartment;
        sP.volume_tt_compartment_rsh = ones(ns,nsh) .* sP.volume_tt_compartment';
        sP.volume_tt_compartment_vec = sP.volume_tt_compartment_rsh(:);
        sP.surface_area_tt_compartment = fitData.surface_area_tt_compartment;
        
        sP.H2O_perm_factor_t = sP.gen_params(5) .* ones(ns,nsh) .* fitData.surface_area_tt_compartment' .* fitData.H2O_molar_volume / 10^9;    % (Liters * cm^3) / (millimoles * ms)
        
    end
    
    sP.num_segs = fitData.uIpms.gen.num_segments;
    sP.num_shells = fitData.uIpms.fit.num_shells;
    
    sP.intra_volume = fitData.intra_volume;
    sP.seg_intra_volume = fitData.seg_intra_volume;
    sP.extra_volume = fitData.extra_volume;
    
    sP.voltage_probe_seg = fitData.uIpms.gen.voltage_probe_segment;
    sP.voltage_probe_shell = fitData.uIpms.gen.voltage_probe_shell;
    sP.current_probe_seg = fitData.uIpms.gen.current_probe_segment;
    sP.current_probe_shell = fitData.uIpms.gen.current_probe_shell;


    sP.Re2 = fitData.uIpms.gen.Re2;    % MOhms
    sP.tau_mu = fitData.uIpms.gen.tau_mu;  % ms
    sP.tau_in = fitData.uIpms.gen.tau_in;    % ms
    sP.mu_gain = fitData.uIpms.gen.mu_gain;  % gain
    sP.R0 = fitData.uIpms.gen.R0;     % output resistance
    
    all_params;      % uncomment so one can watch fitting
    
    
    sP.tdvars = fitData.tdvars;
    sP.tdvars_data = fitData.tdvars_data;
    
    outstruct = [];
    
    for nprot = 1:fitData.num_protocols
        
        clear ss_vals
        
        sP.input_method = fitData.prot{nprot}.input;
        segment_matrices;
        
        if sP.track_volumes
            sP.osmo = fitData.prot{nprot}.conc.osmolarity;     % mM

            sP.fixed_intra_anion_n = fitData.prot{nprot}.conc.fixed_intra_anion_n;   % nanomoles
            sP.fixed_intra_valence = fitData.prot{nprot}.conc.zX;
            sP.neut_osmo_intrac_n = fitData.prot{nprot}.conc.neut_osmo_intrac_n;

            sP.imperm_extra_anion_n = fitData.prot{nprot}.conc.imperm_extra_anion_n; %nanomoles
            sP.neut_osmo_extrac_n = fitData.prot{nprot}.conc.neut_osmo_extrac_n;
            sP.imperm_extra_valence = fitData.prot{nprot}.conc.z_eA;

            if nsh > 0
                sP.imperm_t_anion_n = fitData.prot{nprot}.conc.imperm_t_anion_n';
                sP.neut_osmo_t_n = fitData.prot{nprot}.conc.neut_osmo_t_n';
            end            
        end
        
        sP.Molecs_base = fitData.prot{nprot}.conc.Molecs_base;
        
        if nsh > 0
            for nms = 1:fitData.num_molecs
                if fitData.uIpms.use_old_wrong_diffusion
                    sP.molec_diff_mult{nms} = fitData.prot{nprot}.conc.Molecs_base(nms).diffusion_constant * fitData.uIpms.fit.sigma_t / (10^3 * fitData.delta_r^2);
                else
                    sP.molec_diff_mult{nms} = 2.* fitData.prot{nprot}.conc.Molecs_base(nms).diffusion_constant * fitData.uIpms.fit.sigma_t ./ ( 10^3 * fitData.delta_r .* (sP.radius_tt_for_diff(:,2:end).^2 - sP.radius_tt_for_diff(:,1:end-1).^2) ); % 1/(cm ms)
                end
            end
            sP.use_old_wrong_diffusion = fitData.uIpms.use_old_wrong_diffusion;
        end
        
        if sP.track_icec_ion_numbers
            sP.electrode_cation = fitData.prot{nprot}.electrode_cation;
            sP.electrode_cation_frac = fitData.prot{nprot}.electrode_cation_frac;
            sP.electrode_anion = fitData.prot{nprot}.electrode_anion;
            sP.electrode_anion_frac = fitData.prot{nprot}.electrode_anion_frac;
        end

        if isfield(fitData.prot{nprot},'sig_step_rise')
            sP.sig_step_rise = fitData.prot{nprot}.sig_step_rise;
            sP.sig_step_midtime = fitData.prot{nprot}.sig_step_midtime;
        end
        
        for ncs = 1:fitData.uIpms.num_chans
            sP.chan_partial_cond_factor(ncs) = sP.chan_partial_cond_factor(ncs) * fitData.prot{nprot}.chan_perm_factors(ncs);
            if nsh > 0
                sP.chant_partial_cond_factor{ncs} = sP.chant_partial_cond_factor{ncs} * fitData.prot{nprot}.chan_perm_factors(ncs);
            end
        end
        

        
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % find equilibrium values for beginning of simulation
        
        num_steps = length(fitData.prot{nprot}.steps_to_run);
        
        if isfield(fitData.prot{nprot},'saved_ss_vals') && ~isempty(fitData.prot{nprot}.saved_ss_vals)
            init_vals = fitData.prot{nprot}.saved_ss_vals;
        else

            if fitData.prot{nprot}.uniform_initial_steady_state
                Vm_0 = mean(fitData.prot{nprot}.mean_initial_voltage);  
                Vt_0 = mean(fitData.prot{nprot}.mean_initial_voltage);
                I_ss = mean(fitData.prot{nprot}.mean_initial_current);
                run_init_ss = 1;
            else
                Vm_0 = fitData.prot{nprot}.mean_initial_voltage;  
                Vt_0 = fitData.prot{nprot}.mean_initial_voltage;
                I_ss = fitData.prot{nprot}.mean_initial_current;
                run_init_ss = num_steps;
            end

            init_vals = zeros(fitData.tdvars_data(2),run_init_ss);

            sP.tshift = fitData.prot{nprot}.tshift;

            for j = 1:run_init_ss

                if fitData.uIpms.each_step_diff_exp == 1
                    for ncs = 1:fitData.uIpms.num_chans
                        param_matrix = sP.chan_params_vc{nprot}{ncs};
                        sP.chan_params{ncs}(2:size(param_matrix,1)) = param_matrix(2:end,j);
                    end
                end

                init_vals_ss = [];

                init_vals_ss = [init_vals_ss; Vm_0(j).*ones(ns,1)];
                init_vals_ss = [init_vals_ss; I_ss(j)*fitData.uIpms.gen.R0 + Vm_0(j)];
                init_vals_ss = [init_vals_ss; 0];
                init_vals_ss = [init_vals_ss; Vt_0(j).*ones(ns*nsh,1)];

                if nsh == 0
                    fitData.volume_tt_compartment = 0;
                end
                ttvolumes = ones(ns,nsh).*fitData.volume_tt_compartment';
                init_vals_ss = [init_vals_ss; ttvolumes(:)];
                init_vals_ss = [init_vals_ss; fitData.seg_intra_volume.*ones(ns,1)];
                init_vals_ss = [init_vals_ss; fitData.extra_volume];

                for nms = 1:fitData.num_molecs
                    if fitData.uIpms.track_volumes      % molecule state numbers represents nanomoles (not concentration)
                        if nsh > 0
                            init_tt_moln = reshape(fitData.prot{nprot}.conc.Molecs_base(nms).init_tt .* ones(ns,nsh) .* fitData.volume_tt_compartment' .* 1000,[],1);    
                            init_vals_ss = [init_vals_ss; init_tt_moln];
                        end
                        init_i_moln = fitData.prot{nprot}.conc.Molecs_base(nms).init_intrac .* ones(ns,1) .* fitData.seg_intra_volume .* 1000;
                        init_e_moln = fitData.prot{nprot}.conc.Molecs_base(nms).init_extrac .* fitData.extra_volume .* 1000;
                        init_vals_ss = [init_vals_ss; init_i_moln];
                        init_vals_ss = [init_vals_ss; init_e_moln];
                    else
                        init_vals_ss = [init_vals_ss; fitData.prot{nprot}.conc.Molecs_base(nms).init_tt.*ones(ns*nsh,1)];  % tt values
                        init_vals_ss = [init_vals_ss; fitData.prot{nprot}.conc.Molecs_base(nms).init_intrac.*ones(ns,1)];  % intrac values
                        init_vals_ss = [init_vals_ss; fitData.prot{nprot}.conc.Molecs_base(nms).init_extrac.*ones(1,1)];  % extrac values
                    end

                end

                for ncs = 1:fitData.uIpms.num_chans
                    chan_initvals = fitData.uIpms.Chan_base(ncs).vars{3};
                    for ivc = 1:fitData.uIpms.Chan_base(ncs).vars{1}
                        init_vals_ss = [init_vals_ss; chan_initvals(ivc).*ones(ns,1)]; % for sarc membrane
                        init_vals_ss = [init_vals_ss; chan_initvals(ivc).*ones(ns*nsh,1)]; % for tt membrane
                    end
                end

                sP.times_to_sim = [0;fitData.prot{nprot}.time_to_establish_steady_state];
                sP.times_for_interp = [0;fitData.prot{nprot}.time_to_establish_steady_state];

                if fitData.prot{nprot}.input == 0 || fitData.prot{nprot}.input == 6
                    sP.this_volt_trace_for_interp = Vm_0(j).*[1; 1];
                elseif fitData.prot{nprot}.input == 1 || fitData.prot{nprot}.input == 2 || fitData.prot{nprot}.input == 7
                    sP.this_current_trace_for_interp = I_ss(j).*[1; 1];
                end

                if fitData.prot{nprot}.input == 5 || fitData.prot{nprot}.input == 3
                    sP.Vc_init = fitData.prot{nprot}.sig_step_init(j);
                    sP.Vc_mag = fitData.prot{nprot}.sig_step_final(j);
                else
                    sP.Vc_init = Vm_0(j);
                    sP.Vc_mag = Vm_0(j);
                end

                Vm_all_segs_wparams = @(t,z)Vm_all_segs(t, z, sP);

                [t,u] = ode15s(Vm_all_segs_wparams, sP.times_to_sim, init_vals_ss);

                if sP.track_volumes
                    for ntds = 1:fitData.tdvars_data(1)
                        vec_range = [fitData.tdvars{ntds}(2) : fitData.tdvars{ntds}(3)];
                        ss_vals{ntds}(:,j,:,:) = reshape(u(:, vec_range), [size(u,1), fitData.tdvars{ntds}(4), fitData.tdvars{ntds}(5)]);
                    end
                else
                    for ntds = 1:fitData.tdvars_data(1)
                        vec_range = [fitData.tdvars{ntds}(2) : fitData.tdvars{ntds}(3)];
                        init_vals(vec_range, j) = u(end, vec_range);
                    end
                end

            end

            if sP.track_volumes
                [~, ss_vals, ~, ~] = track_vol_calcs(7, ss_vals, sP, [], nprot);
                for j = 1:run_init_ss
                    for ntds = 1:fitData.tdvars_data(1)
                        vec_range = [fitData.tdvars{ntds}(2) : fitData.tdvars{ntds}(3)];
                        ss_vals_2d = reshape(ss_vals{ntds}(:,j,:,:),[size(u,1), fitData.tdvars{ntds}(4)*fitData.tdvars{ntds}(5)]);
                        init_vals(vec_range,j) = ss_vals_2d(end,:);
                    end
                end

            end

            if 1
                figure(999)
                if sP.track_volumes
                    plot(t,squeeze(ss_vals{1}(:,1,:)))     %Vm
                else
                    plot(t,u(:,1))
                end
            end

        end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % do main simulation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        sP.times_to_sim = fitData.prot{nprot}.times_to_sim;
        num_output_times = length(fitData.prot{nprot}.times_to_sim);

        sP.delta_t = sP.times_to_sim(2) - sP.times_to_sim(1);

        use_ode_time_points = 0;    % 1 only works if only 1 step

        if ~use_ode_time_points
            
            pred_vals = [];
            for ntds = 1:fitData.tdvars_data(1)
                pred_vals{ntds} = zeros(num_output_times, num_steps, fitData.tdvars{ntds}(4), fitData.tdvars{ntds}(5));
            end
            
        end


        for i = 1:num_steps
            
            if fitData.uIpms.each_step_diff_exp == 1
                for ncs = 1:fitData.uIpms.num_chans
                    param_matrix = sP.chan_params_vc{nprot}{ncs};
                    sP.chan_params{ncs}(2:size(param_matrix,1)) = param_matrix(2:end,i);
                end
            end
            

            if fitData.prot{nprot}.input == 0 || fitData.prot{nprot}.input == 6
                sP.times_for_interp = fitData.prot{nprot}.times_for_interp;
                sP.this_volt_trace_for_interp = fitData.prot{nprot}.volt_traces_for_interp(:,i);
            elseif fitData.prot{nprot}.input == 1
                sP.times_for_interp = fitData.prot{nprot}.times_for_interp;
                sP.this_current_trace_for_interp = fitData.prot{nprot}.currents_for_interp(:,i);
            elseif fitData.prot{nprot}.input == 2 || fitData.prot{nprot}.input == 7
                sP.times_for_interp = fitData.prot{nprot}.ic_mat(:,1);
                sP.this_current_trace_for_interp = fitData.prot{nprot}.ic_mat(:,2);
            end

            if fitData.prot{nprot}.input == 1 || fitData.prot{nprot}.input == 2 || fitData.prot{nprot}.input == 7
                if isfield(fitData.prot{nprot},'molecs_change')
                    sP.mols_changing = fitData.prot{nprot}.molecs_change.modified_molecs;
                    for im = sP.mols_changing
                        sP.molch(im).times_for_interp = fitData.prot{nprot}.mc_mat{im}(:,1);
                        sP.molch(im).dmdt_for_interp = fitData.prot{nprot}.mc_mat{im}(:,2);
                    end
                elseif isfield(sP,'molch')
                    sP = rmfield(sP,'molch');
                end
            end
            

            if fitData.prot{nprot}.input == 5
                sP.Vc_mag = fitData.prot{nprot}.sig_step_final(i);
            elseif fitData.prot{nprot}.input == 4
                sP.Vc_mag = fitData.prot{nprot}.voltage_steps(i);
            end

            if fitData.prot{nprot}.uniform_initial_steady_state
                init_vals_use = init_vals(:,1);
            else
                init_vals_use = init_vals(:,i);
            end
            
            optionsode = odeset('RelTol',fitData.uIpms.fit.RelTol,'AbsTol',fitData.uIpms.fit.AbsTol);

            u = init_vals_use';
            tt = sP.times_to_sim(1);

            if fitData.prot{nprot}.input == 2 || fitData.prot{nprot}.input == 7
                Vm_all_segs_wparams = @(t,z)Vm_all_segs(t, z, sP);
                
                time_evals = fitData.prot{nprot}.forced_time_evals;
                num_intervals = length(time_evals) + 1;
                if num_intervals == 1
                    [t,u] = ode15s(Vm_all_segs_wparams, sP.times_to_sim, init_vals_use, optionsode);
                else
                    thistimestart = sP.times_to_sim(1);
                    for nint = 1:num_intervals
                        if nint == num_intervals
                            thistimeend = sP.times_to_sim(end);
                        else
                            thistimeend = time_evals(nint);
                        end
                        simtimes = sP.times_to_sim(sP.times_to_sim >= thistimestart & sP.times_to_sim <= thistimeend);
                        
                        if ismember(thistimeend, sP.times_to_sim)
                            keep_end = 1;
                        else
                            keep_end = 0;
                        end
                        
                        simtimes = union(thistimestart, simtimes);
                        simtimes = simtimes(:);
                        
                        if ~keep_end
                            simtimes = [simtimes; thistimeend];
                        end
                        
                        [~,utmp] = ode15s(Vm_all_segs_wparams, simtimes, init_vals_use, optionsode);
                        
                        init_vals_use = utmp(end,:)';
                        
                        if ~keep_end
                            utmp(end,:) = [];
                        end
                        
                        u = [u; utmp(2:end,:)];
                        
                        thistimestart = thistimeend;
                    end
                end
                tt = sP.times_to_sim;
                    
            else
            
                Vm_all_segs_wparams = @(t,z)Vm_all_segs(t, z, sP);
                [t,u] = ode15s(Vm_all_segs_wparams, sP.times_to_sim, init_vals_use, optionsode);
                tt = sP.times_to_sim;
            end
            
            if ~use_ode_time_points
                num_t = num_output_times;
            else
                num_t = length(tt);
            end
            
            for ntds = 1:fitData.tdvars_data(1)
                vec_range = [fitData.tdvars{ntds}(2) : fitData.tdvars{ntds}(3)];
                try
                    pred_vals{ntds}(:,i,:,:) = reshape(u(:, vec_range), [num_t, fitData.tdvars{ntds}(4), fitData.tdvars{ntds}(5)]);
                catch
                    pred_vals{ntds}(:,i,:,:) = zeros(num_t, fitData.tdvars{ntds}(4), fitData.tdvars{ntds}(5));
                end
                    
            end
        end

        if sP.track_volumes
            [pred_vals, ~, outstruct, state_num] = track_vol_calcs(7, pred_vals, sP, outstruct, nprot);
            
        else
            
            state_num = 7;
            for nms = 1:sP.num_molecs
                state_num = state_num + 1;
                outstruct.prot{nprot}.molecs_t_rsh{nms} = pred_vals{state_num};
                state_num = state_num + 1;
                outstruct.prot{nprot}.molecs_i{nms} = pred_vals{state_num};
                state_num = state_num + 1;
                outstruct.prot{nprot}.molecs_e{nms} = pred_vals{state_num};
            end
            
        end
        
        [Im_curr, It_curr_rsh, ImT, ImC, ItC, Vabs, ImTot, ItTot] = Im_from_Vm(pred_vals, sP, use_ode_time_points, tt);
        
        outstruct.prot{nprot}.ImTot = ImTot;    % sum of current through sarc membrane (cap plus channels) plus t-tubule membrane (cap plus channels) for each segment
        outstruct.prot{nprot}.ItTot = ItTot;
        outstruct.prot{nprot}.Itot_m_pred = sum(ImTot,3);
        outstruct.prot{nprot}.ImC = ImC;
        outstruct.prot{nprot}.Itot_m_C = sum(ImC,3);
        outstruct.prot{nprot}.ImT = ImT;    % sum of current through t-tubule membrane (cap plus channels) for each segment
        outstruct.prot{nprot}.Itot_m_T = sum(ImT,3);
        outstruct.prot{nprot}.ItC = ItC;
        outstruct.prot{nprot}.Itot_t_C = sum(ItC,[3,4]);

        
        for ncs = 1:sP.num_chans
            outstruct.prot{nprot}.ImChan{ncs} = Im_curr{ncs};
            outstruct.prot{nprot}.Itot_m_chan{ncs} = sum(Im_curr{ncs},3);
            if nsh > 0
                outstruct.prot{nprot}.ItChan{ncs} = It_curr_rsh{ncs};
                outstruct.prot{nprot}.Itot_t_chan{ncs} = sum(It_curr_rsh{ncs},[3,4]);
            end
        end
            
        outstruct.prot{nprot}.Vabs = Vabs;
        outstruct.prot{nprot}.Vm_pred = pred_vals{1};
        outstruct.prot{nprot}.V0_pred = pred_vals{2};
        outstruct.prot{nprot}.dV0_dt_pred = pred_vals{3};
        outstruct.prot{nprot}.Vt_pred = pred_vals{4};
        outstruct.prot{nprot}.volume_t_pred = pred_vals{5};
        outstruct.prot{nprot}.volume_i_pred = pred_vals{6};
        outstruct.prot{nprot}.volume_e_pred = pred_vals{7};
        outstruct.prot{nprot}.volume_t_tot = sum(pred_vals{5},[3,4]);
        outstruct.prot{nprot}.volume_i_tot = sum(pred_vals{6},3);
        
        
        for ncs = 1:sP.num_chans
            for nsvs = 1:sP.chan_data(ncs).vars{1}
                state_num = state_num + 1;
                outstruct.prot{nprot}.chan_memb_pred(ncs).sv{nsvs} = pred_vals{state_num};
                state_num = state_num + 1;
                outstruct.prot{nprot}.chan_tt_pred(ncs).sv{nsvs} = pred_vals{state_num};
            end
        end
        
        outstruct.prot{nprot}.simtimes = tt;
        outstruct.prot{nprot}.saved_ss_vals = init_vals;

        if fitData.uIpms.add_capac_coupling
            Vcc = V_from_capac_coupling(tt,fitData.prot{nprot},fitData.uIpms.gen.num_segments,fitData.uIpms.capac_coupling);
            outstruct.prot{nprot}.Vmeas = Vabs + Vcc;
        end
        
    end

end