function [out1, out2] = Chan_Cl_DiFranco(states, V, params, flag, max_cond_factor, use_GHK, RTF, Cl, zCl, varargin)
% states are Clh

%global RTF

out1 = [];
out2 = [];

if flag == 1
    out1 = {'Clh'};
    suggested_init_vals = 0.1;
    out2 = {1, 8, suggested_init_vals, 4, [3]};  % [num_states, num_params, init_vals, num voltage dep outputs, steady-state state output indexes (from flag4,out1)]
    
else
    
    %P_Cl_max = params(1);
    V_alpha_Cl_bar = params(2);
    K_alpha_Cl = params(3); 
    alpha_Cl_bar = params(4);
    V_beta_Cl_bar = params(5);
    K_beta_Cl = params(6);
    beta_Cl_bar = params(7);

    Cl_m = params(8);

    %zCl = -1;

    Vprime_Cl = RTF .* log(Cl_m ./ Cl{1,2}) ./ (-zCl);
    
end

if flag == 3 || flag == 2

    Clh = states{1};
    
    I_Cl_base = find_base_current(V, max_cond_factor, Cl{1,2}, Cl{1,1}, use_GHK, Vprime_Cl, zCl, RTF);
    out2 = I_Cl_base .* Clh;
end
    
if flag == 3 || flag == 4
    
    alpha_Clh = alpha_Cl_bar ./ (1 + exp( -(V - V_alpha_Cl_bar) ./ K_alpha_Cl) );
    beta_Clh = beta_Cl_bar ./ (1 + exp( (V - V_beta_Cl_bar) ./ K_beta_Cl) );

    if flag == 3
        out1{1} = alpha_Clh.*(1-Clh) - beta_Clh.*Clh;
    elseif flag == 4
        out1{1} = alpha_Clh;
        out1{2} = beta_Clh;
        out1{3} = alpha_Clh ./ (alpha_Clh + beta_Clh);
        out1{4} = 1 ./ (alpha_Clh + beta_Clh);
        
        out2 = {'alpha\_Clh','beta\_Clh','Clh\_inf','tau\_Clh'};
    end
    
end
    


