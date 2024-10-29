function [out1, out2] = Chan_NaK_Wallinga(states, V, params, flag, max_cond_factor, use_GHK, RTF, C, z, varargin)
% max_cond_factor is actually a current for this channel
% in C, 1st row is Na, 2nd row is K, 1st col is intrac, 2nd col is extrac

%global RTF

out1 = [];
out2 = [];

if flag == 1
    suggested_init_vals = [];
    out2 = {0, 8, suggested_init_vals, 1, []};  % [num_states, num_params, init_vals vector, num voltage dep outputs]
    
else
    
    %J_NaK= params(1);
    Km_K = params(2);
    Km_Na = params(3);
    c1 = params(4);     % 0.12
    c2 = params(5);     % 0.01
    c3 = params(6);     % 0.04
    c4 = params(7);     % 7
    c5 = params(8);     % 67.3
    
    
    sigma = 1./c4 .* (exp(C{1,2}./c5) - 1);   % 2 3 4
    f = 1./(1 + c1 .* exp(-c2.*V./RTF) + c3 .* sigma .* exp(-V./RTF));  % 2 3 4
    
    if flag == 3 || flag == 2
    
        denom = (1 + Km_K./C{2,2}).^2 .* (1 + Km_Na./C{1,1}).^3; % 2 3
        I_NaK_V_indep = max_cond_factor ./ denom; % 2 3
        out2 = I_NaK_V_indep .* f;                                                 % 2 3
    else
        out1 = {f}; % 4
        out2 = {'f(v)'};    % 4
    end

end


