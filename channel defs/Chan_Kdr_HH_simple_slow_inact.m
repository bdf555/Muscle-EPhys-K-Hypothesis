%function [out1, out2] = Chan_Kdr_HH_simple_slow_inact(states, V, params, flag, max_cond_factor, K_o, K_i, use_GHK, RTF)
function [out1, out2] = Chan_Kdr_HH_simple_slow_inact(states, V, params, flag, max_cond_factor, use_GHK, RTF, K, z, varargin)

out1 = [];
out2 = [];

if flag == 1
    out1 = {'n', 'hK'};
    suggested_init_vals = [0, 1];
    out2 = {2, 10, suggested_init_vals, 6, [3 5]};  % [num_states, num_params, init_vals vector, num voltage dep outputs, num voltage dependent outputs, steady-state state output indexes (from flag4,out1)]
    
else
    
    %P_Kdr_max = params(1);
    V_n_bar = params(2);
    K_alpha_n = params(3);
    K_beta_n = params(4);
    alpha_n_bar = params(5);
    beta_n_bar = params(6);
    
    %S params
    % follows Wallinga 1999
    % dhK/dt = (hK_inf - hK)/tau
    % hK_inf = 1/(1 + exp((V - V50)/Aslope))
    % tau = exp(-(V - tau_num_const) / tau_denom_const)
    V50 = params(7);
    Aslope = params(8);
    tau_num_const = params(9);   
    tau_denom_const = params(10);  
    
end

if flag == 3 || flag == 2

    n = states{1};
    hK = states{2};
    
    I_Kdr_base_memb = find_base_current(V, max_cond_factor, K{1,2}, K{1,1}, use_GHK, 0, z, RTF);
    out2 = I_Kdr_base_memb .* n.^4 .* hK;
end

if flag == 3 || flag == 4
    alpha_n = alpha_n_bar .* (V - V_n_bar) ./ (1 - exp(-(V - V_n_bar)./K_alpha_n));
    beta_n = beta_n_bar .* exp(-(V - V_n_bar)./K_beta_n);
    
    tau = 1000 .* exp(-(V - tau_num_const) ./ tau_denom_const);
    hK_inf = 1./(1 + exp((V - V50)./Aslope));
    
    if flag == 3
        out1{1} = alpha_n.*(1-n) - beta_n.*n;
        out1{2} = (hK_inf - hK) ./ tau;
    elseif flag == 4
        out1 = {alpha_n, beta_n};
        out1{3} = alpha_n ./ (alpha_n + beta_n);
        out1{4} = 1./ (alpha_n + beta_n);
        out1{5} = hK_inf;
        out1{6} = tau;
        
        out2 = {'alpha\_n','beta\_n','n\_inf','tau\_n','hK\_inf','tau\_hK'};
    end
    
end

