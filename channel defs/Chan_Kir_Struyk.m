function [out1, out2] = Chan_Kir_Struyk(states, V, params, flag, max_cond_factor, use_GHK, RTF, C, valence, varargin)

% Based on Struyk and Cannon, 2007, Muscle and Nerve
% Kir V50 shifts with EK

% an HHstate is like the m or h state in Nav in HH formulation
% an 'activating' state goes from closed (0) to open (1) as voltage goes more positive
% an 'inactivating' state goes from open (1) to closed (0) as voltage goes more positive

% for activating:
% S_inf = 1 / (1 + exp(-(V-V50)/Aslope) )

% for inactivating:
% S_inf = 1 / (1 + exp((V-V50)/Aslope) )

% tau is a Gaussian shape over voltage:
% tau = tau_mag * exp(- ((V - tau_Vmean)/tau_Vwidth)^2 )

% rate of approach to S_inf:
% dS/dt = (S_inf - S) / tau

% current:
% I = Ibase * S1^exponent1 * S2^exponent2

out1 = [];
out2 = [];

if flag == 1
    out1 = {'S1', 'S2'};
    suggested_init_vals = [0, 1];   % should reset these according to activating/inactivating type, in parameter file that calls this with flag 1
                                    % example: NaPIC_vars{3} = [0.9, 0];    (where NaPIC_vars is where out2 is assigned)
    out2 = {2, 16, suggested_init_vals, 4, [1 3]};  % [num_states, num_params, init_vals vector, num voltage dep outputs, num voltage dependent outputs, steady-state state output indexes (from flag4,out1)]
    
else
    
    %P_Kdr_max = params(1);
    delta_V50_S1 = params(2);   % When (Vm - EK) equals this delta, 50% of max permeability occurs
    Aslope_S1 = params(3);
    tauS1_mag = params(4);
    tauS1_Vmean = params(5);
    tauS1_Vwidth = params(6);
    
    delta_V50_S2 = params(7);
    Aslope_S2 = params(8);
    tauS2_mag = params(9);
    tauS2_Vmean = params(10);
    tauS2_Vwidth = params(11);
    
    % params below are not meant to be fit
    valence = params(12);   % valence of ion
    S1_act = params(13);    % -1 for activating, 1 for inactivating
    S1_exp = params(14);    % exponent when determining final open fraction (such as the 3 in m^3)
    S2_act = params(15);
    S2_exp = params(16);

end

if flag == 3 || flag == 2

    S1 = states{1};
    S2 = states{2};
    
    I_base_memb = find_base_current(V, max_cond_factor, C{1,2}, C{1,1}, use_GHK, 0, valence, RTF);
    out2 = I_base_memb .* S1.^S1_exp .* S2.^S2_exp;
end

if flag == 3 || flag == 4
    
    Nernst = RTF.*log(C{1,2}./C{1,1})./valence;
    
    %init = -S1_V50shift_maxshift;
    %final = S1_V50shift_maxshift;
    %S1_V50_shift = sigmoid_step_func(Nernst,S1_V50shift_rise,S1_V50shift_offset,final,init);
    %V50_S1_effective = V50_S1 + S1_V50_shift;
    
    S1_inf = 1 ./ (1 + exp(S1_act.*(V-Nernst-delta_V50_S1)./Aslope_S1) );
    tau_S1 = tauS1_mag .* exp(- ((V - tauS1_Vmean)./tauS1_Vwidth).^2 );
    
    %init = -S2_V50shift_maxshift;
    %final = S2_V50shift_maxshift;
    %S2_V50_shift = sigmoid_step_func(Nernst,S2_V50shift_rise,S2_V50shift_offset,final,init);
    %V50_S2_effective = V50_S2 + S2_V50_shift;
    
    S2_inf = 1 ./ (1 + exp(S2_act.*(V-Nernst-delta_V50_S2)./Aslope_S2) );
    tau_S2 = tauS2_mag .* exp(- ((V - tauS2_Vmean)./tauS2_Vwidth).^2 );
    
    if flag == 3
        out1{1} = (S1_inf - S1) ./ tau_S1;
        out1{2} = (S2_inf - S2) ./ tau_S2;
    elseif flag == 4
        out1{1} = S1_inf;
        out1{2} = tau_S1;
        out1{3} = S2_inf;
        out1{4} = tau_S2;
        
        out2 = {'S1\_inf','tau\_S1','S2\_inf','tau\_S2'};
    end
    
end

