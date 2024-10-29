% radial cable model based on Wallinga et. al. (1999) Eur. Biophys. J.
% with corrected K diffusion scaling

function d_dt = Vm_all_segs(t, z, sP)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% map calculated variables to z()

if t > 16959.48
%    t;
    a = 1;    % for breakpoints
end

ns = sP.num_segs;
nsh = sP.num_shells;
vps = sP.voltage_probe_seg;
vpsh = sP.voltage_probe_shell;
ips = sP.current_probe_seg;
ipsh = sP.current_probe_shell;


for ntdv = 1:sP.tdvars_data(1)
    vec_range = [sP.tdvars{ntdv}(2) : sP.tdvars{ntdv}(3)];
    tdv{ntdv} = z(vec_range);
end

Vm = tdv{1};
V0vps = tdv{2};
dV0vps_dt = tdv{3};  % first time deriv of Vm in segment where probes are inserted
Vt = tdv{4};
Vt_rsh = reshape(Vt,[ns,nsh]);

volume_t = tdv{5};
volume_t_rsh = reshape(volume_t,[ns,nsh]);
volume_i = tdv{6};
volume_e = tdv{7};

state_num = 7;

molecs_t = [];
molecs_t_rsh = [];
molecs_i = [];
molecs_e = [];
if ~sP.track_volumes  % states are molecule concentrations
    for nms = 1:sP.num_molecs
        state_num = state_num + 1;
        molecs_t{nms} = tdv{state_num};
        molecs_t_rsh{nms} = reshape(molecs_t{nms},[ns,nsh]);
        state_num = state_num + 1;
        molecs_i{nms} = tdv{state_num};
        state_num = state_num + 1;
        molecs_e{nms} = tdv{state_num};
    end
else  % states are molecule numbers
    molecs_t_n = [];
    molecs_t_rsh_n = [];
    molecs_i_n = [];
    molecs_e_n = [];
    for nms = 1:sP.num_molecs
        state_num = state_num + 1;
        molecs_t_n{nms} = tdv{state_num};
        molecs_t_rsh_n{nms} = reshape(molecs_t_n{nms},[ns,nsh]);
        state_num = state_num + 1;
        molecs_i_n{nms} = tdv{state_num};
        state_num = state_num + 1;
        molecs_e_n{nms} = tdv{state_num};
    end
end    

chan_memb = [];
chan_tt = [];
chan_tt_rsh = [];
for ncs = 1:sP.num_chans
    if sP.chan_data(ncs).vars{1} ~= 0
        for nsvs = 1:sP.chan_data(ncs).vars{1}
            state_num = state_num + 1;
            chan_memb(ncs).sv{nsvs} = tdv{state_num};
            state_num = state_num + 1;
            chan_tt(ncs).sv{nsvs} = tdv{state_num};
            chan_tt_rsh(ncs).sv{nsvs} = reshape(chan_tt(ncs).sv{nsvs},[ns,nsh]);
        end
    else
        chan_memb(ncs).sv = [];
        chan_tt_rsh(ncs).sv = [];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if tracking volume, calculate Vm and Vt from V = Q/C.  Also volume changes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if sP.track_volumes

    [volume_t_rsh, volume_i, volume_e, osmo_t_mobile, osmo_i_mobile, osmo_e_mobile] = Volume_from_Ion_nums(sP, molecs_t_rsh_n, molecs_i_n, molecs_e_n);    
    
    [Vm, Vt_rsh, molecs_t_rsh, molecs_i, molecs_e] = Volt_and_Conc_from_Ion_Nums(sP, molecs_t_rsh_n, molecs_i_n, molecs_e_n, volume_t_rsh, volume_i, volume_e);
    dvolumet_dt = 0;    % cm^3/ms
    dvolumee_dt = 0;
    dvolumei_dt = 0;
else
    osmo_t_mobile = [];
    osmo_i_mobile = [];
    osmo_e_mobile = [];
    dvolumet_dt = 0;    % cm^3/ms
    dvolumee_dt = 0;
    dvolumei_dt = 0;
end    


Im_rates = [];
Im_curr = [];
for ncs = 1:sP.num_chans
    tc = sP.chan_data(ncs);
    mol_array = [];
    for nms = 1:length(tc.molecs_needed)
        molnum = tc.molecs_needed(nms);
        mol_array{nms,1} = molecs_i{molnum};
        mol_array{nms,2} = molecs_e{molnum};
    end
        
    z = sP.Molecs_base(tc.ion_current_mult(1,1)).valence;
    max_cond_factor = find_max_cond_factor(tc.conductance_constant, sP.chan_partial_cond_factor(ncs), mol_array{1,1}, mol_array{1,2}, z, tc.use_GHK, sP.F, sP.RTF);
    [Im_rates{ncs}, Im_curr{ncs}] = tc.ChanFunc(chan_memb(ncs).sv, Vm, sP.chan_params{ncs}, 3, max_cond_factor, tc.use_GHK, sP.RTF, mol_array, z, tc.extra_params);
    
end

if nsh > 0
    
    It_rates = [];
    for ncs = 1:sP.num_chans
        tc = sP.chan_data(ncs);
        mol_array = [];
        for nms = 1:length(tc.molecs_needed)
            molnum = tc.molecs_needed(nms);
            mol_array{nms,1} = molecs_i{molnum};
            mol_array{nms,2} = molecs_t_rsh{molnum};
        end

        z = sP.Molecs_base(tc.ion_current_mult(1,1)).valence;
        max_cond_factor = find_max_cond_factor(tc.conductance_constant, sP.chant_partial_cond_factor{ncs}, mol_array{1,1}, mol_array{1,2}, z, tc.use_GHK, sP.F, sP.RTF);
        [It_rates{ncs}, It_curr_rsh{ncs}] = tc.ChanFunc(chan_tt_rsh(ncs).sv, Vt_rsh, sP.chan_params{ncs}, 3, max_cond_factor, tc.use_GHK, sP.RTF, mol_array, z, tc.extra_params);
    end
    
end

[V, I_input, I_input_mols, dVm_dt_inp, dmol_ex_dt] = find_sim_inputs(t, sP, V0vps, 1, 0, Vm(vps));
% for current clamp, put specified current into probe segment, half in adding +K, half in removing Cl.  
% Balance by doing opposite in extracellular compartment.  Implemented in find_sim_inputs.
% for voltage clamp, determine difference between simmed membrane voltage (Vm) and control voltage (set by input method).  
% Multiply by gain (and imagine dividing by a resistor) to generate voltage clamping current.  Input half as K and half as Cl.  Implemented in find_sim_inputs.


membrane_ionic_current = 0;
for ncr = 1:sP.num_chans
    membrane_ionic_current = membrane_ionic_current + Im_curr{ncr};
end


if sP.input_method == 3
    % rig section
    % todo, fix for IC shells
    % I think only thing that needs changing is set Vvps to specific IC
    % compartment (segment plus shell, not just a segment).  So move this section to after finding all the absolute V
    % values.
    Imvps = membrane_ionic_current(vps) + sP.C*dV0vps_dt;

    Vvps = Vm(vps) + Imvps/sP.g_S;
    [Vc, dVc_dt] = sigmoid_step_func(t, sP.sig_step_rise, sP.sig_step_midtime, sP.Vc_mag, sP.Vc_init);
    part1 = sP.mu_gain*(sP.tau_in*dVc_dt + Vc - Vvps);
    part2 = dV0vps_dt*(sP.tau_in + sP.tau_mu) + V0vps;
    d2V0vps_dt2 = 1./(sP.tau_mu*sP.tau_in).*(part1 - part2);
else
    dV0vps_dt = 0;
    d2V0vps_dt2 = 0;
end


if sP.track_volumes
    dVmdt = 0;
    IT_partial = [];
else
    if ~sP.use_IC_comps

        if nsh > 0
            IT_partial = 1./sP.ra.*(Vm - Vt((nsh-1)*ns+1:ns*nsh) + 1./sP.g_S.*membrane_ionic_current);
        else
            IT_partial = 0;
        end
    end

    if sP.use_IC_comps
        if (ns >= 2 || nsh >= 2)

            d = assemble_d_vector(sP, Vt, Vm, V, I_input);

            Vabs_and_dvtdt = sP.A\d;

            if sP.input_method == 0 || sP.input_method == 4 || sP.input_method == 5
                Vabs = Vabs_and_dvtdt(1:ns*nsh-1);  % Vabs is absolute value of voltage in intracellular compartments
                if ~sP.try_small_matrix
                    dVt_dt = Vabs_and_dvtdt(ns*nsh:2*ns*nsh-1);
                end
                vpssh_idx = (vpsh-1)*ns+vps;
                Vabs = [Vabs(1:vpssh_idx-1); V; Vabs(vpssh_idx:end)];
            else
                Vabs = Vabs_and_dvtdt(1:ns*nsh);  % Vabs is absolute value of voltage in intracellular compartments
                if ~sP.try_small_matrix
                    dVt_dt = Vabs_and_dvtdt(ns*nsh+1:2*ns*nsh);
                end
            end

            dVmdt = 1./sP.C.*( Vabs((nsh-1)*ns+1:ns*nsh).*sP.g_S - Vm.*sP.g_S - membrane_ionic_current);
        else
            if sP.input_method == 0 || sP.input_method == 4 || sP.input_method == 5
                Vabs = V;
            else
                Vabs = Vm + I_input./sP.g_S;
            end

            dVmdt = 1./sP.C.*( Vabs.*sP.g_S - Vm.*sP.g_S - membrane_ionic_current);
        end


    else
        if sP.num_segs >= 2

            d1 = sP.B*Vm;
            d2 = 1./sP.g_S.*sP.B * membrane_ionic_current;
            d = d1 + d2 + membrane_ionic_current + IT_partial;

            if sP.input_method == 0 || sP.input_method == 4 || sP.input_method == 5 || sP.input_method == 6 % voltage trace
                series_voltage_change = sP.g_S./sP.C.*(V - Vm(vps));
                dVmps_dt = series_voltage_change - 1./sP.C.*membrane_ionic_current(vps);

                % remove 1 equation because dVmps_dt is known.
                d(vps-1) = d(vps-1) - sP.d_eqn_corr_above*dVmps_dt;
                d(vps+1) = d(vps+1) - sP.d_eqn_corr_below*dVmps_dt;
                d(vps) = [];
            else % current as input
                d(ips) = d(ips) - I_input;
            end

            dVmdt = sP.A\d;     % all the math occurs here

            if sP.input_method == 0 || sP.input_method == 4 || sP.input_method == 5 || sP.input_method == 6
                dVmdt = [dVmdt(1:vps-1); dVmps_dt; dVmdt(vps:end)];
            end

        else
            if sP.input_method == 0 || sP.input_method == 4 || sP.input_method == 5 || sP.input_method == 6
                dVmdt = 1./sP.C.*(sP.g_S.*(V - Vm) - membrane_ionic_current);
            else
                if nsh > 0
                    div_factor = 1+1./(sP.g_S.*sP.ra);
                else
                    div_factor = 1;
                end
                dVmdt = 1./(div_factor.*sP.C).*(I_input - membrane_ionic_current - IT_partial);
            end

        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%if t>5e7
%    a = 1;    % for breakpoints
%end


if nsh > 0
    if ~sP.use_IC_comps
        [dVt_dt, t_memb_plus_elec_current, tdiffusion_molrate, tmembrane_current, tdiffusion_tt_to_ec_molrate, telectrical_current_flowdown] = tt_rate_calculations_rsh(0, Vt_rsh, It_curr_rsh, sP, Vm, dVmdt, IT_partial, molecs_t_rsh, molecs_e, osmo_t_mobile, osmo_e_mobile);
    else
        if ~sP.try_small_matrix
            [~, t_memb_plus_elec_current, tdiffusion_molrate, tmembrane_current, tdiffusion_tt_to_ec_molrate, telectrical_current_flowdown] = tt_rate_calculations_rsh(0, Vt_rsh, It_curr_rsh, sP, Vm, dVmdt, 0, molecs_t_rsh, molecs_e, osmo_t_mobile, osmo_e_mobile);
        else
            [dVt_dt, t_memb_plus_elec_current, tdiffusion_molrate, tmembrane_current, tdiffusion_tt_to_ec_molrate, telectrical_current_flowdown] = tt_rate_calculations_rsh(Vabs, Vt_rsh, It_curr_rsh, sP, Vm, dVmdt, 0, molecs_t_rsh, molecs_e, osmo_t_mobile, osmo_e_mobile);
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% concentration changes in intracellular and extracellular volumes
% ignore diffusion in the intracellular volume
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_sarcmembrane_current = [];
i_ttmembrane_current = [];
dmoli_dt = [];
dmole_dt = [];
if sP.use_electrical_drift
    diffVm = diff(Vm);  % positive when higher number segments are more positive
    diffidx =  (diffVm>=0) + [1:ns-1]';

    Tconciuse = osmo_i_mobile(diffidx);
end

for nms = 1:sP.num_molecs
    
    if sP.track_icec_ion_numbers

        % concentration changes due to membrane current
        i_sarcmembrane_current{nms} = zeros(ns,1);
        i_sarcmembrane_current{nms}(ips) = I_input_mols{nms};

        for nchs = 1:size(sP.molec_chan{nms},1)
            thischan = sP.molec_chan{nms}(nchs,1);
            thismult = sP.molec_chan{nms}(nchs,2);
            i_sarcmembrane_current{nms} = i_sarcmembrane_current{nms} - Im_curr{thischan}.*thismult;  % nA, from this molecule's movement through channels, negative of normal membrane current definition
        end
        e_sarcmembrane_current_ec_to_ic{nms} = sum(i_sarcmembrane_current{nms});    % nA
        if nsh > 0
            i_ttmembrane_current{nms} = -sum(tmembrane_current{nms}, 2);   % nA
            if sP.use_electrical_drift
                e_electrical_current_ec_to_tt = sum(telectrical_current_flowdown{nms}(:,end));
            end
        else
            i_ttmembrane_current{nms} = 0;
        end

        % concentration changes due to longitudinal electrical drift voltage gradient
        if sP.use_electrical_drift
            moleciuse = molecs_i{nms}(diffidx);
            curr_moleci_flowdown = [sP.gIL .* moleciuse.*diffVm./Tconciuse; 0];
            curr_moleci_flowup = [0; -curr_moleci_flowdown(1:end-1)];

            i_electric_current = curr_moleci_flowdown + curr_moleci_flowup; %nA

        else
        end


        % find molecular rates of change (nanomoles/ms) for each compartment

        dmole_dt{nms} = dmol_ex_dt(nms);
        dmole_dt{nms} = dmole_dt{nms} - e_sarcmembrane_current_ec_to_ic{nms} ./ (sP.F .* sP.Molecs_base(nms).valence .* 10^3);     % nanomoles/ms
        dmoli_dt{nms} = i_sarcmembrane_current{nms}./ (sP.F .* sP.Molecs_base(nms).valence .* 10^3);
        if sP.use_electrical_drift
            dmoli_dt{nms} = dmoli_dt{nms} + i_electric_current ./ (sP.F .* sP.Molecs_base(nms).valence .* 10^3);
        end

        if nsh>0
            dmole_dt{nms} = dmole_dt{nms} + sum(tdiffusion_tt_to_ec_molrate{nms});
            if sP.use_electrical_drift
                dmole_dt{nms} = dmole_dt{nms} - e_electrical_current_ec_to_tt ./ (sP.F .* sP.Molecs_base(nms).valence .* 10^3);
            end
            dmoli_dt{nms} = dmoli_dt{nms} + i_ttmembrane_current{nms} ./ (sP.F .* sP.Molecs_base(nms).valence .* 10^3);
        end

        % if not tracking volumes, convert molecular rates of change to concentration rates of change
        if ~sP.track_volumes
            dmoli_dt{nms} = dmoli_dt{nms} ./ (volume_i .* 10^3);    % mM/ms
            dmole_dt{nms} = dmole_dt{nms} ./ (volume_e .* 10^3);
        end

    else
        dmoli_dt{nms} = zeros(ns,1);
        dmole_dt{nms} = 0;
    end
    
    if sP.track_tt_ion_numbers
        if nsh > 0
            dmolt_dt{nms} = t_memb_plus_elec_current{nms} ./ (sP.F .* sP.Molecs_base(nms).valence .* 10^3) + tdiffusion_molrate{nms};
        else
            dmolt_dt{nms} = 0;
        end
        
        if ~sP.track_volumes
            dmolt_dt{nms} = dmolt_dt{nms} ./ (volume_t_rsh .* 10^3);
        end

    else
        dmolt_dt{nms} = zeros(ns,nsh);
    end
    
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rate equations

d_dt = zeros(sP.tdvars_data(2),1);

d_dt(sP.tdvars{1}(2) : sP.tdvars{1}(3)) = dVmdt;
d_dt(sP.tdvars{2}(2) : sP.tdvars{2}(3)) = dV0vps_dt;
d_dt(sP.tdvars{3}(2) : sP.tdvars{3}(3)) = d2V0vps_dt2;
if nsh > 0
    d_dt(sP.tdvars{4}(2) : sP.tdvars{4}(3)) = dVt_dt(:); % Vt
    d_dt(sP.tdvars{5}(2) : sP.tdvars{5}(3)) = dvolumet_dt(:); % t volume
end
d_dt(sP.tdvars{6}(2) : sP.tdvars{6}(3)) = dvolumei_dt(:); % i volume
d_dt(sP.tdvars{7}(2) : sP.tdvars{7}(3)) = dvolumee_dt(:); % e volume

state_num = 7;
for nms = 1:sP.num_molecs
    state_num = state_num + 1;
    if nsh > 0
        d_dt(sP.tdvars{state_num}(2) : sP.tdvars{state_num}(3)) = dmolt_dt{nms}(:);  % t
    end
    state_num = state_num + 1;
    d_dt(sP.tdvars{state_num}(2) : sP.tdvars{state_num}(3)) = dmoli_dt{nms}(:);  % i
    state_num = state_num + 1;
    d_dt(sP.tdvars{state_num}(2) : sP.tdvars{state_num}(3)) = dmole_dt{nms};  % e
end

for ncs = 1:sP.num_chans
    for nsvs = 1:sP.chan_data(ncs).vars{1}
        state_num = state_num + 1;
        d_dt(sP.tdvars{state_num}(2) : sP.tdvars{state_num}(3)) = Im_rates{ncs}{nsvs}(:);
        state_num = state_num + 1;
        if nsh > 0
            d_dt(sP.tdvars{state_num}(2) : sP.tdvars{state_num}(3)) = It_rates{ncs}{nsvs}(:);
        end
    end
end


