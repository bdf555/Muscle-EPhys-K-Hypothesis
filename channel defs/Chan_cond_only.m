function [out1, out2] = Chan_cond_only(states, V, params, flag, max_cond_factor, use_GHK, RTF, C, z, varargin)
% If use_GHK = 0, then equivalent to constant (non voltage dependent) conductance.  
% if use_GHK = 1, then this is constant permeability (GHK IV curve).  Still no explicit voltage dependence though.  Just the voltage dependence that occurs through GHK equation.

out1 = [];
out2 = [];

if flag == 1

    suggested_init_vals = [];
    out2 = {0, 1, suggested_init_vals, 0, []};  % [num_states, num_params, init_vals vector, num voltage dep outputs, steady-state state output indexes (from flag4,out1)]
    
elseif flag == 3 || flag == 2
    
    %P_Kdr_max = params(1);
    
    I_base_memb = find_base_current(V, max_cond_factor, C{1,2}, C{1,1}, use_GHK, 0, z, RTF);
    out2 = I_base_memb;

end

