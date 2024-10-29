function [out1, out2] = Chan_NaKCl_cotransport_sat(states, V, params, flag, max_cond_factor, use_GHK, RTF, C, z, varargin)
% need to create 3 channels (using this same function), one each for Na and K moving in, and one for Cl moving in
% conductance_constant (in ChanClass) will be 1 for Na and K, and 2 for Cl
% ion_current_mult = [1 1] for Na, [2 1] for the K part, and [3 1] for the Cl part
% the single parameter (like a permeability) will have crazy units of microAmps/(mM^4*cm^2)
% typical value around 3.7e-8?

% enforces a maximum current by limiting the effect of the concentration difference
% to mimic the effect of some maximum number of channels, that is, a saturation effect

%global RTF

out1 = [];
out2 = [];

if flag == 1
    suggested_init_vals = [];
    out2 = {0, 1, suggested_init_vals, 0, []};  % [num_states, num_params, init_vals vector, num voltage dep outputs]
    
else
    
    %Perm like = params(1);
    delta_C_max = params(2);    % the maximum "driving force" due to concentration difference
    delta_C_50 = params(3);
    %c1 = params(4);     % 0.12
    %c2 = params(5);     % 0.01
    %c3 = params(6);     % 0.04
    %c4 = params(7);     % 7
    %c5 = params(8);     % 67.3
    
    
    %sigma = 1./c4 .* (exp(C{1,2}./c5) - 1);   % 2 3 4
    %f = 1./(1 + c1 .* exp(-c2.*V./RTF) + c3 .* sigma .* exp(-V./RTF));  % 2 3 4
    
    if flag == 3 || flag == 2
    
        %denom = (1 + Km_K./C{2,2}).^2 .* (1 + Km_Na./C{1,1}).^3; % 2 3
        %I_NaK_V_indep = max_cond_factor ./ denom; % 2 3
        %out2 = I_NaK_V_indep .* f;                                                 % 2 3
        conc_diff = C{1,2}.*C{2,2}.*C{3,2}.^2 - C{1,1}.*C{2,1}.*C{3,1}.^2;
        
        eff_conc_diff = delta_C_max .* conc_diff ./ (conc_diff + delta_C_50);
        
        %out2 = -max_cond_factor .* conc_diff .* z .* ones(size(V));                                                 % 2 3
        out2 = -max_cond_factor .* eff_conc_diff .* z .* ones(size(V));                                                 % 2 3
    %else
        %out1 = {f}; % 4
        %out2 = {'f(v)'};    % 4
    end

end


