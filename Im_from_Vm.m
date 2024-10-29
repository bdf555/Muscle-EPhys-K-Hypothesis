function [Im_curr, It_curr_rsh, ImT, ImC, ItC, Vabs, ImTot, ItTot] = Im_from_Vm(pred_vals, sP, use_ode_t, tt)
% membrane current from Vm and Vt, voltage across membranes
% assumes 1st dimension of data in Vm, m, and h is time
    
Vm = pred_vals{1};
Vt_all = pred_vals{4};
    
state_num = 7;

molecs_t_rsh = [];
molecs_i = [];
molecs_e = [];
for nms = 1:sP.num_molecs
    state_num = state_num + 1;
    molecs_t_rsh{nms} = pred_vals{state_num};
    state_num = state_num + 1;
    molecs_i{nms} = pred_vals{state_num};
    state_num = state_num + 1;
    molecs_e{nms} = pred_vals{state_num};
end
    

chan_memb = [];
chan_tt_rsh = [];
Im_curr = [];

for ncs = 1:sP.num_chans
    if sP.chan_data(ncs).vars{1} ~= 0
        for nsvs = 1:sP.chan_data(ncs).vars{1}
            state_num = state_num + 1;
            chan_memb(ncs).sv{nsvs} = pred_vals{state_num};
            state_num = state_num + 1;
            chan_tt_rsh(ncs).sv{nsvs} = pred_vals{state_num};
        end
    else
        chan_memb(ncs).sv = [];
        chan_tt_rsh(ncs).sv = [];
    end
    
    tc = sP.chan_data(ncs);
    mol_array = [];
    for nms = 1:length(tc.molecs_needed)
        molnum = tc.molecs_needed(nms);
        mol_array{nms,1} = molecs_i{molnum};
        mol_array{nms,2} = molecs_e{molnum};
    end
        
    z = sP.Molecs_base(tc.ion_current_mult(1,1)).valence;
    max_cond_factor = find_max_cond_factor(tc.conductance_constant, sP.chan_partial_cond_factor(ncs), mol_array{1,1}, mol_array{1,2}, z, tc.use_GHK, sP.F, sP.RTF);
    [~, Im_curr{ncs}] = tc.ChanFunc(chan_memb(ncs).sv, Vm, sP.chan_params{ncs}, 2, max_cond_factor, tc.use_GHK, sP.RTF, mol_array, z, tc.extra_params);

end
    
if ~use_ode_t
    if size(Vm,2) == 1 && ndims(Vm) == 2
        dVmdt = gradient(Vm, sP.delta_t);
    elseif size(Vm,3) == 1
        [~, dVmdt] = gradient(Vm, sP.delta_t);
    else
        [~, dVmdt, ~] = gradient(Vm, sP.delta_t);
    end
else
    if size(Vm,2) == 1 && ndims(Vm) == 2
        dVmdt = gradient(Vm, tt);
    elseif size(Vm,3) == 1
        [~, dVmdt] = gradient(Vm, 1, tt);
    else
        [~, dVmdt, ~] = gradient(Vm, 1, tt, 1);
    end
end    

ImC = sP.C.*dVmdt;


membrane_ionic_current = 0;
for ncr = 1:sP.num_chans
    membrane_ionic_current = membrane_ionic_current + Im_curr{ncr};
end
    
ImTot = ImC + membrane_ionic_current;

ntimes = size(Vt_all,1);
nsteps = size(Vt_all,2);
nsegs = size(Vt_all,3);
nshells = size(Vt_all,4);

ImT = 0;
It_curr_rsh = [];
    
if ~isempty(Vt_all)
    Vt_top = Vt_all(:,:,:,end);
    ImT = 1./sP.ra.*(Vm - Vt_top + 1./sP.g_S.*(sP.C.*dVmdt + membrane_ionic_current));

    ImTot = ImTot + ImT;
    
    if ~use_ode_t
        if size(Vm,2) == 1 && size(Vm,3) == 1 && size(Vt_all,4) == 1
            dVtdt = gradient(Vt_all, sP.delta_t);
        else
            [~, dVtdt] = gradient(Vt_all, sP.delta_t);
        end
    else
        if size(Vm,2) == 1 && size(Vm,3) == 1 && size(Vt_all,4) == 1
            dVtdt = gradient(Vt_all, tt);
        else
            [~, dVtdt] = gradient(Vt_all, 1, tt,1,1);
        end
    end


    for k = 1:size(dVtdt,4)
        ItC(:,:,:,k) = sP.C_t((k-1)*nsegs+1)*dVtdt(:,:,:,k);
    end

    cond_mat = ones(ntimes, nsteps, nsegs, nshells);
    cond_mat = permute(cond_mat,[4 1 2 3]);
    
    for ncs = 1:sP.num_chans
        part_cond = cond_mat .* sP.chant_partial_cond_factor{ncs}(1,:)';
        part_cond = permute(part_cond,[2 3 4 1]);
        
        tc = sP.chan_data(ncs);
        mol_array = [];
        for nms = 1:length(tc.molecs_needed)
            molnum = tc.molecs_needed(nms);
            mol_array{nms,1} = molecs_i{molnum};
            mol_array{nms,2} = molecs_t_rsh{molnum};
        end

        z = sP.Molecs_base(tc.ion_current_mult(1,1)).valence;
        max_cond_factor = find_max_cond_factor(tc.conductance_constant, part_cond, mol_array{1,1}, mol_array{1,2}, z, tc.use_GHK, sP.F, sP.RTF);
        [~, It_curr_rsh{ncs}] = tc.ChanFunc(chan_tt_rsh(ncs).sv, Vt_all, sP.chan_params{ncs}, 2, max_cond_factor, tc.use_GHK, sP.RTF, mol_array, z, tc.extra_params);
    end

    ItTot = ItC;
    for ncs = 1:sP.num_chans
        ItTot = ItTot + It_curr_rsh{ncs};
    end
else
    ItC = 0;
    ItTot = 0;
end

if sP.use_IC_comps

    Vabs = zeros(ntimes,nsteps,nsegs,nshells);

    if sP.input_method == 2

        % complicated because chosen times for when current changes may or
        % may not line up with sP.times_to_sim times, and need to keep
        % results to only include times that match up with data, for
        % fitting

        input_current_table = [sP.times_to_sim(1) sP.input_current_start_value; sP.input_current];
        num_i_vals = size(input_current_table,1);
        for nivals = 1:num_i_vals
            i_input = input_current_table(nivals,2);
            if nivals == num_i_vals
                tend = sP.times_to_sim(end) + 1;
            else
                tend = input_current_table(nivals+1,1);
            end
            t_idxs = find(sP.times_to_sim >= input_current_table(nivals,1) & sP.times_to_sim < tend);
            for nt = t_idxs'
                for nst = 1:nsteps
                    [V, I_input] = find_sim_inputs(sP.times_to_sim(nt), sP, 0, 0, i_input);
                    Vtvec = reshape(Vt_all(nt,nst,:,:),[nsegs*nshells,1]);
                    d = assemble_d_vector(sP, Vtvec, squeeze(Vm(nt,nst,:)), V, I_input);

                    Vabstemp = sP.A(1:nsegs*nshells,1:nsegs*nshells)\d;
                    Vabs(nt,nst,:,:) = reshape(Vabstemp,[nsegs,nshells]);
                end
            end
        end

    else

        for nt = 1:ntimes
            for nst = 1:nsteps
                [V, I_input] = find_sim_inputs(sP.times_to_sim(nt), sP, 0, 0, 0);
                Vtvec = reshape(Vt_all(nt,nst,:,:),[nsegs*nshells,1]);
                d = assemble_d_vector(sP, Vtvec, squeeze(Vm(nt,nst,:)), V, I_input);

                Vabstemp = sP.A(1:nsegs*nshells,1:nsegs*nshells)\d;
                Vabs(nt,nst,:,:) = reshape(Vabstemp,[nsegs,nshells]);
            end
        end
    end
else
    vs = (ImC + membrane_ionic_current)./sP.g_S;
    
    Vabs = Vm + vs;
end

