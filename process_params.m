
uIpms.num_chans = length(Chan_base);

if exist('each_step_diff_exp','var')
    uIpms.each_step_diff_exp = each_step_diff_exp;
else
    uIpms.each_step_diff_exp = 0;
end

if ~exist('current_stim_fit_possible','var')
    gen_params(6) = 1;
    gen_param_bounds = [gen_param_bounds; 1 3 0];
end

if exist('add_capac_coupling','var')
    uIpms.add_capac_coupling = add_capac_coupling;
    if exist('capac_coupling_tpeak','var')
        uIpms.capac_coupling.tpeak = capac_coupling_tpeak;
    else
        uIpms.capac_coupling.tpeak = 0.14;  % ms
    end
    if exist('capac_coupling_doubletau','var')
        uIpms.capac_coupling.doubletau = capac_coupling_doubletau;
    else
        uIpms.capac_coupling.doubletau = 0.2;   % ms
    end
    if exist('capac_coupling_maxmag','var')
        uIpms.capac_coupling.maxmag = capac_coupling_maxmag;
    else
        uIpms.capac_coupling.maxmag = 13;  % mV
    end
    if exist('capac_coupling_falloff','var')
        uIpms.capac_coupling.falloff = capac_coupling_falloff;
    else
        uIpms.capac_coupling.falloff = 70;  % nA
    end
        
else
    uIpms.add_capac_coupling = 0;
end

if exist('use_prev_fit_params') && exist('all_params') && use_prev_fit_params
    uIpms.fit.param_list = all_params;
elseif uIpms.each_step_diff_exp == 0
    uIpms.fit.param_list = [gen_params];
    for nCs = 1:uIpms.num_chans
        uIpms.fit.param_list = [uIpms.fit.param_list Chan_base(nCs).params];
    end
else
    uIpms.fit.param_list = [gen_params];
    paramidxnum = length(gen_params);
    uIpms.fit.param_locator = [];
    for nCs = 1:uIpms.num_chans
        paramidxnum = paramidxnum +1; 
        uIpms.fit.param_list = [uIpms.fit.param_list Chan_base(nCs).params{1}{1}(1)];
        uIpms.fit.param_list_idxnum{nCs}{1}{1}{1} = paramidxnum;
        uIpms.fit.param_locator(paramidxnum,:) = [nCs 1 1 1];
       
        for nprots = 1:length(Chan_base(nCs).params)
            for npns = 2:length(Chan_base(nCs).params{nprots})
                for nstepsnum = 1:length(Chan_base(nCs).params{nprots}{npns})
                    paramidxnum = paramidxnum +1; 
                    uIpms.fit.param_list = [uIpms.fit.param_list Chan_base(nCs).params{nprots}{npns}(nstepsnum)];
                    uIpms.fit.param_list_idxnum{nCs}{nprots}{npns}{nstepsnum} = paramidxnum;
                    uIpms.fit.param_locator(paramidxnum,:) = [nCs nprots npns nstepsnum];
                end
            end
        end
    end
end
    
param_bounds_all = [gen_param_bounds];
if uIpms.each_step_diff_exp == 0
    for nCs = 1:uIpms.num_chans
        param_bounds_all = [param_bounds_all; Chan_base(nCs).param_bounds];
    end
else
    for nCs = 1:uIpms.num_chans
        param_bounds_all_vc{nCs} = Chan_base(nCs).param_bounds;
    end
end

    
thislistnum = length(gen_params);
uIpms.gen_param_idxs = [1 : thislistnum];
cumnum = thislistnum;
uIpms.fit.params_to_fit = [uIpms.gen_param_idxs(fit_gen_params)];

uIpms.molec_chan = [];    % list of which channels affect which molecules, and the multiplier for that molecule
for nms = 1:length(Molecs_base)
    uIpms.molec_chan{nms} = [];
end

if uIpms.each_step_diff_exp == 0
    for nCs = 1:uIpms.num_chans

        thislistnum = length(Chan_base(nCs).params);
        uIpms.chan_param_idxs{nCs} = [cumnum+1 : cumnum + thislistnum];
        cumnum = cumnum + thislistnum;

        for molecn = 1:size(Chan_base(nCs).ion_current_mult,1)
            thismolec = Chan_base(nCs).ion_current_mult(molecn,1);
            uIpms.molec_chan{thismolec} = [uIpms.molec_chan{thismolec}; nCs Chan_base(nCs).ion_current_mult(molecn,2)];
        end

        uIpms.fit.params_to_fit = [uIpms.fit.params_to_fit uIpms.chan_param_idxs{nCs}(Chan_base(nCs).fit_params)];
    end
else
    for nCs = 1:uIpms.num_chans
        for molecn = 1:size(Chan_base(nCs).ion_current_mult,1)
            thismolec = Chan_base(nCs).ion_current_mult(molecn,1);
            uIpms.molec_chan{thismolec} = [uIpms.molec_chan{thismolec}; nCs Chan_base(nCs).ion_current_mult(molecn,2)];
        end
        if Chan_base(nCs).fit_params{1}{1}(1) == 1
            uIpms.fit.params_to_fit = [uIpms.fit.params_to_fit uIpms.fit.param_list_idxnum{nCs}{1}{1}{1}];
        end
        for nprots = 1:length(Chan_base(nCs).params)
            for npns = 2:length(Chan_base(nCs).params{nprots})
                for nstepsnum = 1:length(Chan_base(nCs).params{nprots}{npns})
                    if Chan_base(nCs).fit_params{nprots}{npns}(nstepsnum) == 1
                        uIpms.fit.params_to_fit = [uIpms.fit.params_to_fit uIpms.fit.param_list_idxnum{nCs}{nprots}{npns}{nstepsnum}];
                    end
                end
            end
        end
    end
end
    

if uIpms.each_step_diff_exp == 0
    j = 1;
    uIpms.fit.lb_params = [];
    uIpms.fit.ub_params = [];
    for i = uIpms.fit.params_to_fit
        if param_bounds_all(i,1) == 1
            uIpms.fit.lb_params(j) = uIpms.fit.param_list(i) ./ param_bounds_all(i,2);
            uIpms.fit.ub_params(j) = uIpms.fit.param_list(i) .* param_bounds_all(i,2);
        elseif param_bounds_all(i,1) == 2
            uIpms.fit.lb_params(j) = uIpms.fit.param_list(i) - param_bounds_all(i,2);
            uIpms.fit.ub_params(j) = uIpms.fit.param_list(i) + param_bounds_all(i,2);
        elseif param_bounds_all(i,1) == 3
            uIpms.fit.lb_params(j) = param_bounds_all(i,2);
            uIpms.fit.ub_params(j) = param_bounds_all(i,3);
        end
        j = j+1;
    end
else
    j = 1;
    uIpms.fit.lb_params = [];
    uIpms.fit.ub_params = [];
    for i = uIpms.fit.params_to_fit
        param_val = uIpms.fit.param_list(i);
        if i <= 5
            param_bounds = param_bounds_all(i,:);
        else
            param_data = uIpms.fit.param_locator(i,:);
            param_bounds = param_bounds_all_vc{param_data(1)}(param_data(3),:);
        end
        if param_bounds(1) == 1
            uIpms.fit.lb_params(j) = param_val ./ param_bounds(2);
            uIpms.fit.ub_params(j) = param_val .* param_bounds(2);
        elseif param_bounds(1) == 2
            uIpms.fit.lb_params(j) = param_val - param_bounds(2);
            uIpms.fit.ub_params(j) = param_val + param_bounds(2);
        elseif param_bounds(1) == 3
            uIpms.fit.lb_params(j) = param_bounds(2);
            uIpms.fit.ub_params(j) = param_bounds(3);
        end
        j = j+1;
    end

end

uIpms.gen.Temp = Temperature;

uIpms.Molecs_base = Molecs_base;
uIpms.Chan_base = Chan_base;

uIpms.track_volumes = track_volumes;
uIpms.use_electrical_drift = use_electrical_drift;
if track_volumes
    uIpms.osmolarity = osmolarity;
    uIpms.gen.zX = zX;    % valence of impermeant intracellular anion
    uIpms.gen.z_eA = z_eA;  % valence of eA
end

if track_volumes == 1
    uIpms.track_icec_ion_numbers = 1;
elseif exist('track_icec_ion_numbers','var')
    uIpms.track_icec_ion_numbers = track_icec_ion_numbers;
else
    uIpms.track_icec_ion_numbers = 1;
end

if track_volumes == 1
    uIpms.track_tt_ion_numbers = 1;
elseif exist('track_tt_ion_numbers','var')
    uIpms.track_tt_ion_numbers = track_tt_ion_numbers;
else
    uIpms.track_tt_ion_numbers = 1;
end

if exist('use_old_wrong_diffusion','var')
    uIpms.use_old_wrong_diffusion = use_old_wrong_diffusion;
else
    uIpms.use_old_wrong_diffusion = 0;
end

if exist('chan_ids_for_plot_and_fit','var')
    uIpms.chan_ids_for_plot_and_fit = chan_ids_for_plot_and_fit;
else
    uIpms.chan_ids_for_plot_and_fit = 1;
end
    

uIpms.gen.num_segments = num_segments;
uIpms.gen.voltage_probe_segment = voltage_probe_segment;
uIpms.gen.voltage_probe_shell = voltage_probe_shell;
uIpms.gen.current_probe_segment = current_probe_segment;
uIpms.gen.current_probe_shell = current_probe_shell;
uIpms.gen.cell_radius = cell_radius;
uIpms.gen.cell_length = cell_length;
uIpms.gen.extrac_volume_multiple = extrac_volume_multiple;

uIpms.protocols = protocol;
if exist('filename','var')
    uIpms.filenames = filename;
end
uIpms.fit.do_fits_or_onerun = do_fits_or_onerun;
uIpms.fit.use_IC_comps = use_IC_comps;

uIpms.fit.time_resolution_increase = time_resolution_increase;

uIpms.fit.num_shells = num_shells;
uIpms.fit.zeta = zeta;
uIpms.fit.rho_t = rho_t;
uIpms.fit.sigma_t = sigma_t;

uIpms.fit.use_access_resistance = use_access_resistance;
uIpms.fit.Ra = Ra;

uIpms.fit.fit_algorithm = fit_algorithm;
uIpms.fit.MaxFunctionEvaluations = MaxFunctionEvaluations;
if exist('use_parallel','var')
    uIpms.fit.use_parallel = use_parallel;
else
    uIpms.fit.use_parallel = false;
end
if fit_algorithm == 1
    uIpms.fit.MaxIterations = MaxIterations;
    uIpms.fit.FiniteDifferenceStepSize = FiniteDifferenceStepSize;
    uIpms.fit.FiniteDifferenceType = FiniteDifferenceType;
    uIpms.fit.StepTolerance = StepTolerance;
    uIpms.fit.FunctionTolerance = FunctionTolerance;
elseif fit_algorithm == 3
    uIpms.fit.MinSurrogatePts = MinSurrogatePts;
elseif fit_algorithm == 4
    uIpms.fit.MaxTime = MaxTime;
elseif fit_algorithm == 5
    uIpms.fit.SwarmSize = SwarmSize;
    uIpms.fit.MaxTime = MaxTime;
end

uIpms.fit.RelTol = RelTol;
uIpms.fit.AbsTol = AbsTol;


uIpms.gen.Re2 = Re2;    % MOhms
uIpms.gen.tau_mu = tau_mu;  % ms
uIpms.gen.tau_in = Re1*(C_in - C_n);    % ms
uIpms.gen.mu_gain = mu_gain;  % gain
uIpms.gen.R0 = R0;

uIpms.pfile = pfile;
uIpms.pgfile = pgfile;




