function d = assemble_d_vector(sP, Vtvec, Vmvec, V, I_input)

d = sP.B*Vtvec;

d((sP.num_shells-1)*sP.num_segs+1:sP.num_segs*sP.num_shells) = d((sP.num_shells-1)*sP.num_segs+1:sP.num_segs*sP.num_shells) - Vmvec.*sP.g_S; 
if sP.input_method == 0 || sP.input_method == 4 || sP.input_method == 5 || sP.input_method == 6 % voltage trace
    d = d - sP.d_eqn_corr*V;
    d((sP.voltage_probe_shell-1)*sP.num_segs + sP.voltage_probe_seg) = [];
else
    d((sP.current_probe_shell-1)*sP.num_segs + sP.current_probe_seg) = d((sP.current_probe_shell-1)*sP.num_segs + sP.current_probe_seg) - I_input;
end
