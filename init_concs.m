function protocol = init_concs(protocol_in, fitData)

protocol = protocol_in;
protocol.conc.Molecs_base = fitData.uIpms.Molecs_base;

if isfield(protocol_in, 'electrode_cation')
    cation_num = protocol_in.electrode_cation;
else
    cation_num = 2;
end
if isfield(protocol_in, 'electrode_cation_frac')
    cation_frac = protocol_in.electrode_cation_frac;
else
    cation_frac = 0.5;
end
if isfield(protocol_in, 'electrode_anion')
    anion_num = protocol_in.electrode_anion;
else
    anion_num = 3;
end
anion_frac = 1 - cation_frac;

protocol.electrode_cation = cation_num;
protocol.electrode_cation_frac = cation_frac;
protocol.electrode_anion = anion_num;
protocol.electrode_anion_frac = anion_frac;

if fitData.uIpms.track_volumes
    num_molecs = length(fitData.uIpms.Molecs_base);

    for nm = 1:num_molecs
        Mol_o(nm) = fitData.uIpms.Molecs_base(nm).init_extrac;
        Mol_i(nm) = fitData.uIpms.Molecs_base(nm).init_intrac;
    end

    osmo = fitData.uIpms.osmolarity;
    zX = fitData.uIpms.gen.zX;    % valence of impermeant intracellular anion
    z_eA = fitData.uIpms.gen.z_eA;      % valence of eA

    if isfield(protocol_in, 'conc')
        for nm = 1:num_molecs
            if isfield(protocol_in.conc, 'modified_molecs')
                if ismember(nm,protocol_in.conc.modified_molecs)
                    if isfield(protocol_in.conc.Molecs_base(nm), 'init_extrac') && ~isempty(protocol_in.conc.Molecs_base(nm).init_extrac)
                        Mol_o(nm) = protocol_in.conc.Molecs_base(nm).init_extrac;
                    end
                    if isfield(protocol_in.conc.Molecs_base(nm), 'init_intrac') && ~isempty(protocol_in.conc.Molecs_base(nm).init_intrac)
                        Mol_i(nm) = protocol_in.conc.Molecs_base(nm).init_intrac;
                    end
                end
            end
        end
        if isfield(protocol_in.conc, 'osmolarity')
            osmo = protocol_in.conc.osmolarity;
        end
        if isfield(protocol_in.conc, 'zX')
            zX = protocol_in.conc.zX;
        end
        if isfield(protocol_in.conc, 'z_eA')
            z_eA = protocol_in.conc.z_eA;
        end
    end

    extra_net_charge = 0;
    for nm = 1:num_molecs
        extra_net_charge = extra_net_charge + Mol_o(nm) * fitData.uIpms.Molecs_base(nm).valence;
    end
    %extra_net_charge = Na_o + K_o - Cl_o;     % extra net anion needed
    eA = -extra_net_charge/z_eA;    % nM
    if eA < 0
        error('Extracellular impermeant anion concentration set to negative value due to electroneutrality in format_fit_data_res_inc')
    end
    disp(['Set extracellular impermeant anion to ' num2str(eA) ' mM'])

    intra_net_charge = 0;
    for nm = 1:num_molecs
        intra_net_charge = intra_net_charge + Mol_i(nm) * fitData.uIpms.Molecs_base(nm).valence;
    end
    %intra_net_charge = Na_i + K_i - Cl_i;    % intra net anion needed
    X_i = -intra_net_charge/zX; % nM
    disp(['Set intracellular fixed anion to ' num2str(X_i) ' mM'])

    % use intracellular osmolarity as "true" and fixed value

    neut_osmo_extrac = osmo;
    for nm = 1:num_molecs
        neut_osmo_extrac = neut_osmo_extrac - Mol_o(nm);
    end
    neut_osmo_extrac = neut_osmo_extrac - eA;

    neut_osmo_intrac = osmo;
    for nm = 1:num_molecs
        neut_osmo_intrac = neut_osmo_intrac - Mol_i(nm);
    end
    neut_osmo_intrac = neut_osmo_intrac - X_i;

    %neut_osmo_extrac = osmo - (Na_o + K_o + Cl_o + eA);
    %neut_osmo_intrac = osmo - (Na_i + K_i + Cl_i + X_i);

    if neut_osmo_extrac < 0
        error(['sum of osmolarity of extracellular components is greater than specified osmolarity ('  num2str(osmo) ') in format_fit_data_res_inc.'])
    elseif neut_osmo_intrac < 0
        error(['sum of osmolarity of intracellular components is greater than specified osmolarity ('  num2str(osmo) ') in format_fit_data_res_inc '])
    end

    disp(['Added neutral extracellular osmolyte: ' num2str(neut_osmo_extrac) ' mOsm.']);
    disp(['Added neutral intracellular osmolyte: ' num2str(neut_osmo_intrac) ' mOsm.']);

    C = fitData.uIpms.fit.param_list(1).*fitData.tot_surf_area; %uF



    Q_surplus_intra = protocol.mean_initial_voltage .* C;   %nanoCoulombs
    num_Q_surplus = Q_surplus_intra / fitData.F; % nanomoles

    % make half of surplus Cl, and half (deficit) K
    Mol_i(protocol.electrode_anion) = Mol_i(protocol.electrode_anion) - num_Q_surplus*protocol.electrode_anion_frac/fitData.intra_volume/1e3;
    Mol_i(protocol.electrode_cation) = Mol_i(protocol.electrode_cation) + num_Q_surplus*protocol.electrode_cation_frac/fitData.intra_volume/1e3;

    if fitData.uIpms.fit.num_shells > 0
        for nm = 1:num_molecs
            Mol_t(nm) = Mol_o(nm);
        end

        sarc_area_frac = fitData.sarc_surf_area ./ fitData.tot_surf_area;
        % put ~1/5 of surplus into extracellular space:
        Mol_o(protocol.electrode_anion) = Mol_o(protocol.electrode_anion) + num_Q_surplus*sarc_area_frac*protocol.electrode_anion_frac/fitData.extra_volume/1e3;
        Mol_o(protocol.electrode_cation) = Mol_o(protocol.electrode_cation) - num_Q_surplus*sarc_area_frac*protocol.electrode_cation_frac/fitData.extra_volume/1e3;
        % put ~4/5 of surplus into t-tubule space:
        Mol_t(protocol.electrode_anion) = Mol_t(protocol.electrode_anion) + num_Q_surplus*(1-sarc_area_frac)*protocol.electrode_anion_frac/fitData.tt_total_volume/1e3;
        Mol_t(protocol.electrode_cation) = Mol_t(protocol.electrode_cation) - num_Q_surplus*(1-sarc_area_frac)*protocol.electrode_cation_frac/fitData.tt_total_volume/1e3;
    else
        Mol_o(protocol.electrode_anion) = Mol_o(protocol.electrode_anion) + num_Q_surplus*protocol.electrode_anion_frac/fitData.extra_volume/1e3;
        Mol_o(protocol.electrode_cation) = Mol_o(protocol.electrode_cation) - num_Q_surplus*protocol.electrode_cation_frac/fitData.extra_volume/1e3;
    end

    for nm = 1:num_molecs
        protocol.conc.Molecs_base(nm).init_extrac = Mol_o(nm);
        protocol.conc.Molecs_base(nm).init_intrac = Mol_i(nm);
        if fitData.uIpms.fit.num_shells > 0
            protocol.conc.Molecs_base(nm).init_tt = Mol_t(nm);
        end
    end

    protocol.conc.osmolarity = osmo;
    protocol.conc.fixed_intra_anion_n = X_i .* fitData.seg_intra_volume .* 1000;    % nanomoles of fixed anion in each segment
    protocol.conc.neut_osmo_intrac_n = neut_osmo_intrac .* fitData.seg_intra_volume .* 1000;    % nanomoles of fixed neutral osmolyte in each segment
    if fitData.uIpms.fit.num_shells > 0
        protocol.conc.imperm_t_anion_n = eA .* fitData.volume_tt_compartment .* 1000;  % nanomoles
        protocol.conc.neut_osmo_t_n = neut_osmo_extrac .* fitData.volume_tt_compartment .* 1000;  % nanomoles
    end
    protocol.conc.X_i = X_i;
    protocol.conc.zX = zX;
    protocol.conc.eA = eA;
    protocol.conc.z_eA = z_eA;
    protocol.conc.imperm_extra_anion_n = eA * fitData.extra_volume * 1000;
    protocol.conc.neut_osmo_extrac_n = neut_osmo_extrac * fitData.extra_volume * 1000;
end


