function [V, I_input, I_input_mols, dVm_dt, dmol_ex_dt] = find_sim_inputs(t, sP, V0vps, running_sim, i_input, vmvps)

V = 0;
I_input = 0;
dVm_dt = 0;
for nms = 1:sP.num_molecs
    dmol_ex_dt(nms) = 0;
end

if sP.input_method == 0 || sP.input_method == 6
    V = interp1(sP.times_for_interp, sP.this_volt_trace_for_interp, t - sP.tshift);
    Vslope = gradient(sP.this_volt_trace_for_interp, sP.times_for_interp);
    dVm_dt = interp1(sP.times_for_interp, Vslope, t - sP.tshift);
elseif sP.input_method == 1 || sP.input_method == 2 || sP.input_method == 7
    I_input = interp1(sP.times_for_interp, sP.this_current_trace_for_interp, t - sP.tshift);
    if isfield(sP,'molch')
        for im = sP.mols_changing
            dmol_ex_dt(im) = interp1(sP.molch(im).times_for_interp, sP.molch(im).dmdt_for_interp, t - sP.tshift);
        end
    end

elseif sP.input_method == 3
    I_input = V0vps/sP.R0;
elseif sP.input_method == 4 || sP.input_method == 5
    [V, ~] = sigmoid_step_func(t, sP.sig_step_rise, sP.sig_step_midtime, sP.Vc_mag, sP.Vc_init);
end

% for voltage clamping while tracking volumes:
if sP.track_volumes && sP.track_icec_ion_numbers && (sP.input_method == 0 || sP.input_method == 6 || sP.input_method == 4 || sP.input_method == 5)
    delta_V = V - vmvps;
    I_input = sP.mu_gain*delta_V/sP.R0;
end

for nms = 1:sP.num_molecs
    I_input_mols{nms} = 0;
end

if sP.track_icec_ion_numbers
    % make assumption that half of input current is K ions being pushed in, and half is Cl ions being pulled out (still neg current using convention)
    % a 'pos' stim current is pushing pos ions in through electrode, which is a neg current using conventions, but due to the push of positive ions in, pos ions must now flow out through membrane to keep charge balance...
    I_input_mols{sP.electrode_cation} = I_input*sP.electrode_cation_frac;
    I_input_mols{sP.electrode_anion} = I_input*sP.electrode_anion_frac;
end

