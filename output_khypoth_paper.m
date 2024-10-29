
% let output_mep run first (happens naturally with running of main.m)

% run once for each "type" in 1st section.  Leave all other sections to not run.  This stores important data from simulation.  So this file run 4 times.
% After each type of simulation data stored, run with all of 1st section not running, just to save variables (2nd section)
% Then comment out saving of data and uncomment loading of data and run a final time to make plots.

do_low = 0; % else do high K and hyperKPP

% allow specific section below depending on which condition/mutation was run
if 0
    wt_Vm = outstruct.prot{1}.Vm_pred(:,1,fitData.uIpms.gen.voltage_probe_segment);
    wt_molrate = Molrate_current;
    wt_ImChan = outstruct.prot{1}.ImChan;
    wt_molecsi = outstruct.prot{1}.molecs_i;
    wt_molecse = outstruct.prot{1}.molecs_e;
    wt_all_params = all_params;
elseif 0
    hypo_Vm = outstruct.prot{1}.Vm_pred(:,1,fitData.uIpms.gen.voltage_probe_segment);
    hypo_molrate = Molrate_current;
    hypo_ImChan = outstruct.prot{1}.ImChan;
    hypo_molecsi = outstruct.prot{1}.molecs_i;
    hypo_molecse = outstruct.prot{1}.molecs_e;
    hypo_all_params = all_params;
elseif 0
    hyper_Vm = outstruct.prot{1}.Vm_pred(:,1,fitData.uIpms.gen.voltage_probe_segment);
    hyper_molrate = Molrate_current;
    hyper_ImChan = outstruct.prot{1}.ImChan;
    hyper_molecsi = outstruct.prot{1}.molecs_i;
    hyper_molecse = outstruct.prot{1}.molecs_e;
    hyper_all_params = all_params;
elseif 0
    at_Vm = outstruct.prot{1}.Vm_pred(:,1,fitData.uIpms.gen.voltage_probe_segment);
    at_molrate = Molrate_current;
    at_ImChan = outstruct.prot{1}.ImChan;
    at_molecsi = outstruct.prot{1}.molecs_i;
    at_molecse = outstruct.prot{1}.molecs_e;
    at_all_params = all_params;
end

%save("khypothvars_low.mat","wt_Vm","wt_molrate","wt_ImChan","wt_molecsi","wt_molecse","wt_all_params","hypo_Vm","hypo_molrate","hypo_ImChan","hypo_molecsi","hypo_molecse","hypo_all_params",...
%    "hyper_Vm","hyper_molrate","hyper_ImChan","hyper_molecsi","hyper_molecse","hyper_all_params","at_Vm","at_molrate","at_ImChan","at_molecsi","at_molecse","at_all_params")
%save("khypothvars_high.mat","wt_Vm","wt_molrate","wt_ImChan","wt_molecsi","wt_molecse","wt_all_params","hypo_Vm","hypo_molrate","hypo_ImChan","hypo_molecsi","hypo_molecse","hypo_all_params",...
%    "hyper_Vm","hyper_molrate","hyper_ImChan","hyper_molecsi","hyper_molecse","hyper_all_params","at_Vm","at_molrate","at_ImChan","at_molecsi","at_molecse","at_all_params")

%load("khypothvars_low.mat")
%load("khypothvars_high.mat")


% make plots
if 0
ttimes = fitData.prot{nprot}.times_to_sim/60000 - 500;

all_vars = ttimes;

if do_low
    x_limits = [-10 170];
else
    x_limits = [-20 500];
end

% Vm plot
figure(1010)
clf
plot(ttimes, wt_Vm)
hold on
if do_low
    plot(ttimes, hypo_Vm)
else
    plot(ttimes, hyper_Vm)
end
plot(ttimes, at_Vm)
title('Membrane potential')
if do_low
    legend('wild type, lowK', 'hypoKPP, lowK', 'AT, lowK' )
else
    legend('wild type, highK', 'hyperKPP, highK', 'AT, highK' )
end
xlabel('minutes');
ylabel('mV')
    
xlim(x_limits)

all_vars = [all_vars,wt_Vm,hypo_Vm,hyper_Vm,at_Vm];


% net Na current plot
figure(1011)
clf
plot(ttimes, wt_molrate{1}/total_cap)
hold on
if do_low
    plot(ttimes, hypo_molrate{1}/total_cap)
else
    plot(ttimes, hyper_molrate{1}/total_cap)
end
plot(ttimes, at_molrate{1}/total_cap)
title('Net Na current')
if do_low
    legend('wild type, lowK', 'hypoKPP, lowK', 'AT, lowK' )
else
    legend('wild type, highK', 'hyperKPP, highK', 'AT, highK' )
end
xlim(x_limits)
xlabel('minutes');
ylabel('nA/uF')

all_vars = [all_vars,wt_molrate{1}/total_cap,hypo_molrate{1}/total_cap,hyper_molrate{1}/total_cap,at_molrate{1}/total_cap];

% net K current plot
figure(1012)
clf
plot(ttimes, wt_molrate{2}/total_cap)
hold on
if do_low
    plot(ttimes, hypo_molrate{2}/total_cap)
else
    plot(ttimes, hyper_molrate{2}/total_cap)
end    
plot(ttimes, at_molrate{2}/total_cap)
title('Net K current')
if do_low
    legend('wild type, lowK', 'hypoKPP, lowK', 'AT, lowK' )
else
    legend('wild type, highK', 'hyperKPP, highK', 'AT, highK' )
end
xlim(x_limits)
xlabel('minutes');
ylabel('nA/uF')

all_vars = [all_vars,wt_molrate{2}/total_cap,hypo_molrate{2}/total_cap,hyper_molrate{2}/total_cap,at_molrate{2}/total_cap];

% Kex conc
figure(1013)
clf
plot(ttimes, wt_molecse{2})
hold on
if do_low
    plot(ttimes, hypo_molecse{2})
else
    plot(ttimes, hyper_molecse{2})
end
plot(ttimes, at_molecse{2})
title('Extracellular [K]')
if do_low
    legend('wild type, lowK', 'hypoKPP, lowK', 'AT, lowK' )
else
    legend('wild type, highK', 'hyperKPP, highK', 'AT, highK' )
end
xlim(x_limits)
xlabel('minutes');
ylabel('mM')

all_vars = [all_vars,wt_molecse{2},hypo_molecse{2},hyper_molecse{2},at_molecse{2}];


% Na_int conc
figure(1014)
clf
plot(ttimes, wt_molecsi{1})
hold on
if do_low
    plot(ttimes, hypo_molecsi{1})
else
    plot(ttimes, hyper_molecsi{1})
end    
plot(ttimes, at_molecsi{1})
title('Intracellular [Na]')
if do_low
    legend('wild type, lowK', 'hypoKPP, lowK', 'AT, lowK' )
else
    legend('wild type, highK', 'hyperKPP, highK', 'AT, highK' )
end
xlim(x_limits)
xlabel('minutes');
ylabel('mM')

if do_low
    time_index = 30600;
    legend_string = '10 min after';
else
    time_index = 36000;
    legend_string = '100 min after';
end

all_vars = [all_vars,wt_molecsi{1},hypo_molecsi{1},hyper_molecsi{1},at_molecsi{1}];
    
% Kir vs V, wt
figure(1015)
clf
hold on
% for Kir vs V plots
chosen_time_indices = [29999 time_index];
cs = 2;     % selects Kir channel
ctinum = 0;
vidx = find(V == 0,1);
V(vidx) = 0.01;
all_vvars = V';
for cti = chosen_time_indices

    ctinum = ctinum + 1;

    if ctinum == 1
        colorst = 'k';
    elseif ctinum == 2
        colorst = 'r';
    end
    
    tc = uIpms.Chan_base(cs);

    num_v_outs = tc.vars{4};
    mol_array = [];
    for nms = 1:length(tc.molecs_needed)
        molnum = tc.molecs_needed(nms);
        mol_array{nms,1} = wt_molecsi{molnum}(cti,plotstep,plotseg);
        mol_array{nms,2} = wt_molecse{molnum}(cti,plotstep);
    end
    z = Molecs_base(tc.ion_current_mult(1,1)).valence;
    
    [vouts_chan{cs}.vo, vout_names] = tc.ChanFunc([],V,wt_all_params(uIpms.chan_param_idxs{cs}),4,[],[],fitData.RTF, mol_array, z, tc.extra_params); 

    ss_states_idxs = tc.vars{5};
    chan_states = [];
    for ssnum = 1:length(ss_states_idxs)
        chan_states = [chan_states ,vouts_chan{cs}.vo(ss_states_idxs(ssnum))];
    end
    thischan_params = wt_all_params(uIpms.chan_param_idxs{cs});
    thischan_perm = thischan_params(1);
    thischan_partial_cond_factor = thischan_perm .* fitData.sarc_seg_area .* 1000;
        
    thischan_max_cond_factor = find_max_cond_factor(tc.conductance_constant, thischan_partial_cond_factor, mol_array{1,1}, mol_array{1,2}, z, tc.use_GHK, fitData.F, fitData.RTF);
    [~, ImoChan{cs}] = tc.ChanFunc(chan_states, V, wt_all_params(uIpms.chan_param_idxs{cs}), 2, thischan_max_cond_factor, tc.use_GHK, fitData.RTF, mol_array, z, tc.extra_params);

    plot(V,ImoChan{cs}/total_cap, colorst)

    all_vvars = [all_vvars,ImoChan{cs}'/total_cap];

    xline(wt_Vm(cti,1,fitData.uIpms.gen.voltage_probe_segment), colorst)
    wt_Vm(cti,1,fitData.uIpms.gen.voltage_probe_segment)

end
tempzero = zeros(1,length(V));
plot(V,tempzero,'k')
title('Wild Type Kir vs V')
legend('','baseline','',legend_string,'FontSize',16)
xlabel('mV');
ylabel('nA/uF')



% Kir vs V, hypo
figure(1016)
clf
hold on
% for Kir vs V plots
chosen_time_indices = [29999 time_index];
cs = 2;     % selects Kir channel
ctinum = 0;
vidx = find(V == 0,1);
V(vidx) = 0.01;
for cti = chosen_time_indices

    ctinum = ctinum + 1;

    if ctinum == 1
        colorst = 'k';
    elseif ctinum == 2
        colorst = 'r';
    end
    
    tc = uIpms.Chan_base(cs);

    num_v_outs = tc.vars{4};
    mol_array = [];
    for nms = 1:length(tc.molecs_needed)
        molnum = tc.molecs_needed(nms);
        mol_array{nms,1} = hypo_molecsi{molnum}(cti,plotstep,plotseg);
        mol_array{nms,2} = hypo_molecse{molnum}(cti,plotstep);
    end
    z = Molecs_base(tc.ion_current_mult(1,1)).valence;
    
    [vouts_chan{cs}.vo, vout_names] = tc.ChanFunc([],V,hypo_all_params(uIpms.chan_param_idxs{cs}),4,[],[],fitData.RTF, mol_array, z, tc.extra_params); 

    ss_states_idxs = tc.vars{5};
    chan_states = [];
    for ssnum = 1:length(ss_states_idxs)
        chan_states = [chan_states ,vouts_chan{cs}.vo(ss_states_idxs(ssnum))];
    end
    thischan_params = hypo_all_params(uIpms.chan_param_idxs{cs});
    thischan_perm = thischan_params(1);
    thischan_partial_cond_factor = thischan_perm .* fitData.sarc_seg_area .* 1000;
        
    thischan_max_cond_factor = find_max_cond_factor(tc.conductance_constant, thischan_partial_cond_factor, mol_array{1,1}, mol_array{1,2}, z, tc.use_GHK, fitData.F, fitData.RTF);
    [~, ImoChan{cs}] = tc.ChanFunc(chan_states, V, hypo_all_params(uIpms.chan_param_idxs{cs}), 2, thischan_max_cond_factor, tc.use_GHK, fitData.RTF, mol_array, z, tc.extra_params);

    plot(V,ImoChan{cs}/total_cap, colorst)

    all_vvars = [all_vvars,ImoChan{cs}'/total_cap];

    xline(hypo_Vm(cti,1,fitData.uIpms.gen.voltage_probe_segment), colorst)
    hypo_Vm(cti,1,fitData.uIpms.gen.voltage_probe_segment)

end
tempzero = zeros(1,length(V));
plot(V,tempzero,'k')
title('HypoKPP Kir vs V')
legend('','baseline','',legend_string,'FontSize',16)
xlabel('mV');
ylabel('nA/uF')




% Kir vs V, hyper
figure(1017)
clf
hold on
% for Kir vs V plots
chosen_time_indices = [29999 time_index];
cs = 2;     % selects Kir channel
ctinum = 0;
vidx = find(V == 0,1);
V(vidx) = 0.01;
for cti = chosen_time_indices

    ctinum = ctinum + 1;

    if ctinum == 1
        colorst = 'k';
    elseif ctinum == 2
        colorst = 'r';
    end
    
    tc = uIpms.Chan_base(cs);

    num_v_outs = tc.vars{4};
    mol_array = [];
    for nms = 1:length(tc.molecs_needed)
        molnum = tc.molecs_needed(nms);
        mol_array{nms,1} = hyper_molecsi{molnum}(cti,plotstep,plotseg);
        mol_array{nms,2} = hyper_molecse{molnum}(cti,plotstep);
    end
    z = Molecs_base(tc.ion_current_mult(1,1)).valence;
    
    [vouts_chan{cs}.vo, vout_names] = tc.ChanFunc([],V,hyper_all_params(uIpms.chan_param_idxs{cs}),4,[],[],fitData.RTF, mol_array, z, tc.extra_params); 

    ss_states_idxs = tc.vars{5};
    chan_states = [];
    for ssnum = 1:length(ss_states_idxs)
        chan_states = [chan_states ,vouts_chan{cs}.vo(ss_states_idxs(ssnum))];
    end
    thischan_params = hyper_all_params(uIpms.chan_param_idxs{cs});
    thischan_perm = thischan_params(1);
    thischan_partial_cond_factor = thischan_perm .* fitData.sarc_seg_area .* 1000;
        
    thischan_max_cond_factor = find_max_cond_factor(tc.conductance_constant, thischan_partial_cond_factor, mol_array{1,1}, mol_array{1,2}, z, tc.use_GHK, fitData.F, fitData.RTF);
    [~, ImoChan{cs}] = tc.ChanFunc(chan_states, V, hyper_all_params(uIpms.chan_param_idxs{cs}), 2, thischan_max_cond_factor, tc.use_GHK, fitData.RTF, mol_array, z, tc.extra_params);

    plot(V,ImoChan{cs}/total_cap, colorst)

    all_vvars = [all_vvars,ImoChan{cs}'/total_cap];

    xline(hyper_Vm(cti,1,fitData.uIpms.gen.voltage_probe_segment), colorst)
    hyper_Vm(cti,1,fitData.uIpms.gen.voltage_probe_segment)

end
tempzero = zeros(1,length(V));
plot(V,tempzero,'k')
title('HyperKPP Kir vs V')
legend('','baseline','',legend_string,'FontSize',16)
xlabel('mV');
ylabel('nA/uF')



% Kir vs V, at
figure(1018)
clf
hold on
% for Kir vs V plots
chosen_time_indices = [29999 time_index];
cs = 2;     % selects Kir channel
ctinum = 0;
vidx = find(V == 0,1);
V(vidx) = 0.01;
for cti = chosen_time_indices

    ctinum = ctinum + 1;

    if ctinum == 1
        colorst = 'k';
    elseif ctinum == 2
        colorst = 'r';
    end
    
    tc = uIpms.Chan_base(cs);

    num_v_outs = tc.vars{4};
    mol_array = [];
    for nms = 1:length(tc.molecs_needed)
        molnum = tc.molecs_needed(nms);
        mol_array{nms,1} = at_molecsi{molnum}(cti,plotstep,plotseg);
        mol_array{nms,2} = at_molecse{molnum}(cti,plotstep);
    end
    z = Molecs_base(tc.ion_current_mult(1,1)).valence;
    
    [vouts_chan{cs}.vo, vout_names] = tc.ChanFunc([],V,at_all_params(uIpms.chan_param_idxs{cs}),4,[],[],fitData.RTF, mol_array, z, tc.extra_params); 

    ss_states_idxs = tc.vars{5};
    chan_states = [];
    for ssnum = 1:length(ss_states_idxs)
        chan_states = [chan_states ,vouts_chan{cs}.vo(ss_states_idxs(ssnum))];
    end
    thischan_params = at_all_params(uIpms.chan_param_idxs{cs});
    thischan_perm = thischan_params(1);
    thischan_partial_cond_factor = thischan_perm .* fitData.sarc_seg_area .* 1000;
        
    thischan_max_cond_factor = find_max_cond_factor(tc.conductance_constant, thischan_partial_cond_factor, mol_array{1,1}, mol_array{1,2}, z, tc.use_GHK, fitData.F, fitData.RTF);
    [~, ImoChan{cs}] = tc.ChanFunc(chan_states, V, at_all_params(uIpms.chan_param_idxs{cs}), 2, thischan_max_cond_factor, tc.use_GHK, fitData.RTF, mol_array, z, tc.extra_params);

    plot(V,ImoChan{cs}/total_cap, colorst)

    all_vvars = [all_vvars,ImoChan{cs}'/total_cap];

    xline(at_Vm(cti,1,fitData.uIpms.gen.voltage_probe_segment), colorst)
    at_Vm(cti,1,fitData.uIpms.gen.voltage_probe_segment)

end
tempzero = zeros(1,length(V));
plot(V,tempzero,'k')
title('AT Kir vs V')
legend('','baseline','',legend_string,'FontSize',16)
xlabel('mV');
ylabel('nA/uF')

figure(1019)
clf
plot(ttimes, wt_ImChan{2} / total_cap)
hold on
if do_low
    plot(ttimes, hypo_ImChan{2} / total_cap)
else
    plot(ttimes, hyper_ImChan{2} / total_cap)
end
plot(ttimes, at_ImChan{2} / total_cap)
title('Kir current')
if do_low
    legend('wild type, lowK', 'hypoKPP, lowK', 'AT, lowK' )
else
    legend('wild type, highK', 'hyperKPP, highK', 'AT, highK' )
end
xlabel('minutes');
ylabel('nA/uF')
    
xlim(x_limits)

all_vars = [all_vars,wt_ImChan{2} / total_cap,hypo_ImChan{2} / total_cap,hyper_ImChan{2} / total_cap,at_ImChan{2} / total_cap];


figure(1020)
clf
plot(ttimes, wt_ImChan{4} / total_cap)
hold on
if do_low
    plot(ttimes, hypo_ImChan{4} / total_cap)
else
    plot(ttimes, hyper_ImChan{4} / total_cap)
end
plot(ttimes, at_ImChan{4} / total_cap)
title('NaK current')
if do_low
    legend('wild type, lowK', 'hypoKPP, lowK', 'AT, lowK' )
else
    legend('wild type, highK', 'hyperKPP, highK', 'AT, highK' )
end
xlabel('minutes');
ylabel('nA/uF')
    
xlim(x_limits)
all_vars = [all_vars,wt_ImChan{4} / total_cap,hypo_ImChan{4} / total_cap,hyper_ImChan{4} / total_cap,at_ImChan{4} / total_cap];


figure(1021)
clf
plot(ttimes, wt_ImChan{5} / total_cap)
hold on
if do_low
    plot(ttimes, hypo_ImChan{5} / total_cap)
else
    plot(ttimes, hyper_ImChan{5} / total_cap)
end
plot(ttimes, at_ImChan{5} / total_cap)
title('NaP current')
if do_low
    legend('wild type, lowK', 'hypoKPP, lowK', 'AT, lowK' )
else
    legend('wild type, highK', 'hyperKPP, highK', 'AT, highK' )
end
xlabel('minutes');
ylabel('nA/uF')
    
xlim(x_limits)
all_vars = [all_vars,wt_ImChan{5} / total_cap,hypo_ImChan{5} / total_cap,hyper_ImChan{5} / total_cap,at_ImChan{5} / total_cap];

end    
