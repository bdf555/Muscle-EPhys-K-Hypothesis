function [Vm, Vt, molecs_t_rsh, molecs_i, molecs_e] = Volt_and_Conc_from_Ion_Nums(sP, molecs_t_rsh_n, molecs_i_n, molecs_e_n, volume_t_rsh, volume_i, volume_e)

% voltage in mV, ion numbers in nanomoles, concentrations in millimolar, volumes in microliters, capacitance in nanoFarads
% all variables associated with t-tubules must be 2 (or more)-d matrices with num_segments in 1st matrix dimension and num_shells in second
% variables associated with intracellular compartments must have num_segments in 1st matrix dimension

wv_Q_i = sP.fixed_intra_anion_n .* sP.fixed_intra_valence;  %nanomoles, weighted by valence, converted to net charge below
for nms = 1:sP.num_molecs
    wv_Q_i = wv_Q_i + molecs_i_n{nms}.*sP.Molecs_base(nms).valence;
    molecs_i{nms} = molecs_i_n{nms}./volume_i;
    molecs_e{nms} = molecs_e_n{nms}./volume_e;
end
Q_i = wv_Q_i .* sP.F;   % nC

if sP.num_shells>0
    wv_Q_t = sP.imperm_t_anion_n .* sP.imperm_extra_valence; %nanomoles, weighted by valence, converted to net charge below

    for nms = 1:sP.num_molecs
        wv_Q_t = wv_Q_t + molecs_t_rsh_n{nms}.*sP.Molecs_base(nms).valence;
        molecs_t_rsh{nms} = molecs_t_rsh_n{nms}./volume_t_rsh;
    end
    Q_t = wv_Q_t .* sP.F;   % nC
    Qt_summed_shells = sum(Q_t,2);
    if ndims(Q_t) == 2  % no list of times present, evaluating at a single time
        Qi_sl = Q_i + squeeze(Qt_summed_shells);    % net charge (in nanoCoulombs) in intracellular compartment that is left over for sarcolemma after subtracting excess intracellular charge used on tt membranes
    else    % evaluating at multiple times
        sz_Qtss = size(Qt_summed_shells);
        Qi_sl = Q_i + reshape(Qt_summed_shells, sz_Qtss(1), sz_Qtss(3));    % net charge (in nanoCoulombs) in intracellular compartment that is left over for sarcolemma after subtracting excess intracellular charge used on tt membranes
    end
    Vt = -Q_t .* 1000 ./ sP.C_t;  % in mV
else
    Qi_sl = Q_i;
    Vt = [];
    molecs_t_rsh = [];
end

Vm = Qi_sl .* 1000 ./ sP.C;

