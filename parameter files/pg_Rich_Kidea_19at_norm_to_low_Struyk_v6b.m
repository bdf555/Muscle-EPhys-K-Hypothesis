
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                            
% cable values

if 0
num_segments = 19;   % 1 segment works, 3 or more work.  I haven't checked case of 2 segments.
num_shells = 9;				% number of radial shells for both TT and intracellular
voltage_probe_segment = 2;  % segment where V electrode is placed.  Must be less than or equal to num_segments.
voltage_probe_shell = 5;  % shell where V electrode is placed.  Must be less than or equal to num_shells. (1 is center)
current_probe_segment = 2;  % segment where I electrode is placed.  Must be less than or equal to num_segments.  Not used if voltage input is used (input method 0).
current_probe_shell = 5;  % shell where I electrode is placed.  Must be less than or equal to num_shells.  Not used if voltage input is used (input method 0).
elseif 0
num_segments = 5;   % 1 segment works, 3 or more work.  I haven't checked case of 2 segments.
num_shells = 3;				% number of radial shells for both TT and intracellular
voltage_probe_segment = 2;  % segment where V electrode is placed.  Must be less than or equal to num_segments.
voltage_probe_shell = 2;  % shell where V electrode is placed.  Must be less than or equal to num_shells. (1 is center)
current_probe_segment = 2;  % segment where I electrode is placed.  Must be less than or equal to num_segments.  Not used if voltage input is used (input method 0).
current_probe_shell = 2;  % shell where I electrode is placed.  Must be less than or equal to num_shells.  Not used if voltage input is used (input method 0).
else
num_segments = 1;   % 1 segment works, 3 or more work.  I haven't checked case of 2 segments.
num_shells = 0;				% number of radial shells for both TT and intracellular
voltage_probe_segment = 1;  % segment where V electrode is placed.  Must be less than or equal to num_segments.
voltage_probe_shell = 2;  % shell where V electrode is placed.  Must be less than or equal to num_shells. (1 is center)
current_probe_segment = 1;  % segment where I electrode is placed.  Must be less than or equal to num_segments.  Not used if voltage input is used (input method 0).
current_probe_shell = 2;  % shell where I electrode is placed.  Must be less than or equal to num_shells.  Not used if voltage input is used (input method 0).
end

cell_radius = 30e-4;  % cm
%cell_radius = 20e-4;  % cm
cell_length = 1; % cm
%cell_length = .1; % cm


variable_intracellular_volume = 1;     % intracellular volume is allowed to change to keep osmolarity constant. Was called 'track_volumes' before.
infinite_extracellular_volume = 0;  % if 1, then renders extrac_volume_multiple meaningless.  This value is not used in ode version (track time course)
                                    % since the code has no provision for infinite extracellular volume.  This value is used in the steady state solver version.
track_volumes = variable_intracellular_volume;

extrac_volume_multiple = 0.5;     % extracellular volume = intracellular volume times this multiple
                                % interstitial muscle volume is about 25% (0.25) of the total muscle volume (Kandarian, 1991, J Appl Physiol
%track_volumes = 0;  % for now, must have Na as 1st molecule in Molecs_base, K as 2nd, and Cl as 3rd to use this
use_electrical_drift = 0;   % when currents flow along tube, they are carried by ions.  So ion concentrations can change this way.
                            % Must use if t-tubules are on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% output focus
%plotseg = 1;
plotseg = voltage_probe_segment;
plotshell = voltage_probe_shell;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Temperature = 305;			% K

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
use_RM_orig = 0;        % This is not used for steady state version, since must have initial RM_orig guess.
RM_orig = -60;          % an initial guess for RM (if variable_intracellular_volume) or the starting point of RM
                        % If use_RM_orig = 0, then initial guess for time course version is set to Nernst for K based on initial values

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The values below are only used in the steady state version of the code
fix_RM = 0;
use_multiple_init_guesses = 1;
RM_guess = [-55.7 -70 -85 -100];
Nai_guess = [2 4 8 20.67];
Ki_guess = [60 84.7 100 140];
Cli_guess = [6 10 14.4 20];

if track_volumes
    osmolarity = 300;   % mOsm
    %X_i = 108.5;      % mM intracellular impermeant negative ions
    zX = -1.6477;    % valence of impermeant intracellular anion
    %Z_t = 0.05;        % mM t-tubule fixed anion.  Normal t-tubule Cl concentration will be assumed to be equal to Cl_o plus zZt*Z_t.
    %zZt = -1;       % valence of fixed t-tubule anion
    %eA = 40;    % extracellular impermeant anion concentration (lactate, bicarbonate, proteins) mM.  This is present in all extracellular space (tt plus extra)
                % code assumes that the number (nanomoles) of these anions in each compartment (tt and extra) never changes.  This probably isn't strictly true
                % if the relative volumes between tt and extra space changes, but such a possibility is hard to code.
    z_eA = -1;  % valence of eA
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% T-TUBULE PARAMETERS
%zeta = 1e-6;			% cm, volume-to-surface ratio of tts
zeta = 2.87020e-6;
%rho_t = 0.003;			% fraction of fiber volume occupied by t-tubules
rho_t = 0.0117;
sigma_t = 0.34;			% tortuosity factor: fraction of radial directed tubular branches 

use_access_resistance = 1;
Ra = 0.150;				% kohm-cm2, access resistance per cm2 of surface membrane, 
                        % if Ra not used, access resistance will be set at G_T (gen_params(4))
                        % adjusted for geometry of last tt shell and half of shell-to-shell distance


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


volt = rho_t * pi * cell_length * cell_radius.^2;   % cm^3
voli = pi * cell_length * cell_radius^2 - volt;

%voli = 2.7944e-5;
vole = voli*extrac_volume_multiple;
%vole = 5.5887e-5;
%volt = 3.3081e-7;


if 1    % no ion shifting
Molecs_base(1) = MolecClass;    % Na
Molecs_base(1).init_extrac = 137;    % normal
%Molecs_base(1).init_extrac = 128.5;    % high
%Molecs_base(1).init_extrac = 126.5;    % high   v2
%Molecs_base(1).init_extrac = 125.5;    % high   v3
%Molecs_base(1).init_extrac = 139.5;     % low

%Molecs_base(1).init_intrac = 7;
Molecs_base(1).init_intrac = 13.2;

Molecs_base(1).init_tt = 149;
Molecs_base(1).name = 'Na';
end

if 1    % no ion shifting
Molecs_base(2) = MolecClass;    % K
Molecs_base(2).init_extrac = 13;     % normal
%Molecs_base(2).init_extrac = 21.5;     % high
%Molecs_base(2).init_extrac = 23.5;     % high   v2
%Molecs_base(2).init_extrac = 24.5;     % high   v3
%Molecs_base(2).init_extrac = 10.5;      % low

%Molecs_base(2).init_intrac = 110;
Molecs_base(2).init_intrac = 94.1;

Molecs_base(2).init_tt = 1;
Molecs_base(2).name = 'K';
end    
    
if 1    % no ion shifting
Molecs_base(3) = MolecClass;    % Cl
Molecs_base(3).init_extrac = 110;
Molecs_base(3).init_intrac = 6.5;
Molecs_base(3).init_tt = 110;
Molecs_base(3).name = 'Cl';
Molecs_base(3).valence = -1;
elseif 0    % no ion shifting
Molecs_base(3) = MolecClass;    % Cl
Molecs_base(3).init_extrac = 130;
Molecs_base(3).init_intrac = 6.5;
Molecs_base(3).init_tt = 110;
Molecs_base(3).name = 'Cl';
Molecs_base(3).valence = -1;
else
Molecs_base(3) = MolecClass;    % Cl
Molecs_base(3).init_extrac = conc_Clo_final;
Molecs_base(3).init_intrac = conc_Cli_final;
Molecs_base(3).init_tt = conc_Clo_final;
Molecs_base(3).name = 'Cl';
Molecs_base(3).valence = -1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fittable params


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% general parameters

fit_gen_params = [];

% Parameter 1, uF/cm^2, membrane capacitance per area
                    % area is surface (if no TT) or total (if TT)
%Cm = 0.47e-3;   
gen_params(1)= 0.9;   % specific capacitance, uF/cm^2
                                

% Parameter 2, mS/cm^2, series conductance
%G_S = 100;
%gen_params(2) = 100;    % series conductance, ms/cm^2
gen_params(2) = 100000000;    % series conductance, ms/cm^2
                                
% Parameter 6, kOhm-cm, intracellular medium resistivity
                        % Resistance between 1 seg and next in kOhms = 
                        % rho_i.* (cell_length ./ num_segments) ./ (pi.*cell_radius.^2)
%rho_i = 0.180;
% Parameter 3, mS/cm, intracellular medium conductivity
%G_I = 5.5555556;
%G_I = 50;
%G_I = 10000;
%gen_params(3) = 5.55556;    % intracellular medium conductivity, mS/cm
gen_params(3) = 18;    % intracellular medium conductivity, mS/cm
                            
% Parameter 4, conductivity t-tubule lumen, mS/cm
% 100 Ohm-cm = 10 mS/cm; 200 Ohm-cm = 5 mS/cm
%G_T = 3.7;
gen_params(4) = 8;    % t-tubule lumen conductivity, mS/cm

% parameter 5, water permeability through membrane
gen_params(5) = 0.0156;    % cm/s

% for bounds list below, each row is for entry in parameter list above
% 1st value meaning: 1 if 2nd value is multiplier, 2 if 2nd value is
% additive. For example, [1 3] means lower bound = guess_val/3 and upper bound =
% guess_val*3. [2 3] means lower bound = guess_val-3 and upper bound =
% guess_val+3.
gen_param_bounds = [2 0.4 0;  %1
                1 4 0;    %2
                1 1.3 0;  %3
                1 1.3 0];  %4



chan_number = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nav params
if 0
chan_number = chan_number + 1;
Nav = ChanClass;   % Nav
Nav.name = 'Nav';
Nav.ChanFunc = @Chan_Na_HH_simple_slow_inact;
Nav.molecs_needed = [1];    % according to Molecs_base above.  Also list in same order that ChanFunc expects.
Nav.ion_current_mult = [1 1];   % 1st column is molecule id number, 2nd column is multiplier for that ion current vs net current through this channel.
Nav.eta = 1;    % density of channel in t-tubules relative to its density in sarcolemma

Nav.use_GHK = 1;

Nav.fit_params = [];

%G_Na_max = 115.0;   % mS/cm^2
%G_Na_max = 230.0;   % mS/cm^2
%G_Na_max = 430.0;   % mS/cm^2
%G_Na_max = 0.0;   % mS/cm^2
Na_o = Molecs_base(1).init_extrac;
Na_i = Molecs_base(1).init_intrac;
%P_Na_max = permeability_from_conductance(G_Na_max,Na_o,Na_i,1,Temperature);    % cm/s
P_Na_max = 2e-3;
%Nav_params(1) = 7e-4;   %Nav permeability, cm/s
Nav.params(1) = P_Na_max;   %Nav permeability, cm/s

% Na activation
Nav.params(2) = -39.1;  % V_m_bar, mV
Nav.params(3) = 9.8599; %  K_alpha_m, mV
Nav.params(4) = 34.859;  %  K_beta_m, mV
Nav.params(5) = 0.12984;    % alpha_m_bar, 1/(ms*mV)
Nav.params(6) = 0.75173;       % beta_m_bar, 1/ms

% Na inactivation
Nav.params(7) = -35.596;  %  V_h_bar, mV
Nav.params(8) = 25.071; %  K_alpha_h, mV
Nav.params(9) = 3.7965;  %  K_beta_h, mV
Nav.params(10) = 0.0034315;    %  alpha_h_bar, 1/ms
Nav.params(11) = 8.689;       %  beta_h_bar1, 1/ms
%Nav_params(12) = 0.08;     %  beta_h_bar2, 1/(ms*mV)

Nav.params(12) = -87; % V50, mV
Nav.params(13) = 9.8; % Aslope, mV
Nav.params(14) = 60000; % tau_num, ms
Nav.params(15) = 1.5; % tau_denom_const1, dimensionless
Nav.params(16) = 6; % tau_denom_const2
Nav.params(17) = -90; % tau_denom_Vmid, mV
Nav.params(18) = 70; % tau_denom_const3, mV
Nav.params(19) = 2; % tau_denom_exp, dimensionless

Nav.param_bounds = [1 3 0;    %1
    
                2 30 0;   %2
                1 4 0;    %3
                1 4 0;    %4
                1 4 0;    %5
                1 4 0;    %6
                
                2 30 0;   %7
                1 4 0;    %8
                1 4 0;    %9
                1 4 0;    %10    
                1 4 0;    %11
                
                2 30 0;   %12
                2 5 0;    %13
                1 4 0;    %14
                1 4 0;    %15
                1 4 0;    %16
                2 40 0;   %17
                1 4 0;    %18
                1 2 0];   %19
            
            
%[Nav_statenames, Nav_vars] = Nav_func([],[],[],1,[],[],[],[],[]); 
[Nav.statenames, Nav.vars] = Nav.ChanFunc([],[],[],1,[],[],[],[],[]); 
Chan_base(chan_number) = Nav;

% Nav_names contains names of state variables.  = {'m', 'h', 'S'}
% Nav_vars{1} contains number of state variables.  2 in this case
% Nav_vars{2} contains number of parameters. 12 in this case
% Nav_vars{3} contains vector of initial values for sim. (Note, these are initial values to start the 
%  finding of equilibrium, not the actual sim.)  [0,1] in this case.
% Nav_vars{4} contains the number of voltage dependent outputs used for
%  final plotting
% Nav_vars{5} contains the indexes for steady-state state variables from
%  list of voltage dependent outputs

 
% If one wants to assign initial values for the state variables to something other than the default, 
% assign them as in the following statement for the m and h time dependent variables:
% Nav_vars{3} = [0, 0.5] 

% for bounds lists, each row is for entry in parameter list above
% 1st value meaning: 1 if 2nd value in row is multiplier, 2 if 2nd value in row is
% additive, 3 if 2nd and 3rd values in row are lower and upper bound respectively. 
% For example, [1 3] means lower bound = guess_val/3 and upper bound =
% guess_val*3. [2 3] means lower bound = guess_val-3 and upper bound =
% guess_val+3.  [3 20 50] means lower bound is 20, upper bound is 50.

end                                                %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kdr params
if 1
chan_number = chan_number + 1;
Kdr = ChanClass;   % Kdr
Kdr.name = 'Kdr';
Kdr.ChanFunc = @Chan_Kdr_HH_simple_slow_inact;
Kdr.molecs_needed = [2];
Kdr.ion_current_mult = [2 1];
Kdr.eta = 1;

Kdr.use_GHK = 1;

Kdr.fit_params = [];

%G_KDr_max = 18.5;   % mS/cm^2
%K_o = Molecs_base(2).init_extrac;
%K_i = Molecs_base(2).init_intrac;
%P_KDr_max = permeability_from_conductance(G_KDr_max,K_o,K_i,1,Temperature);    % cm/s
P_KDr_max = 2.864e-04;
Kdr.params(1) = P_KDr_max;   %Kdr permeability, cm/s

%Kdr.params(2) = -43.854;  % V_n_bar, mV
Kdr.params(2) = -49.001;  % V_n_bar, mV
Kdr.params(3) = 4.192; %  K_alpha_n, mV
Kdr.params(4) = 17.42;  %  K_beta_n, mV
Kdr.params(5) = 0.013621;    % alpha_n_bar, 1/(ms*mV)
Kdr.params(6) = 0.083997;       % beta_n_bar, 1/ms

% slow inact
Kdr.params(7) = -40;    % V50, mV
%Kdr_params(7) = 100;    % V50, mV
Kdr.params(8) = 7.5;    % Aslope, mV
Kdr.params(9) = -40;     % tau_num_const, mV
Kdr.params(10) = 25.75; % tau_denom_const, mV

Kdr.param_bounds = [1 3 0;
    
                2 30 0;   
                1 4 0;    
                1 4 0;    
                1 4 0;    
                1 4 0;
                
                2 30 0;
                2 5 0;
                2 30 0;
                1 4 0];   
            
[Kdr.statenames, Kdr.vars] = Kdr.ChanFunc([],[],[],1,[],[],[],[],[]); 
Kdr.vars{3} = [1.9e-5, 0.998]; 
Chan_base(chan_number) = Kdr;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Kir params

if 1
chan_number = chan_number + 1;
Kir = ChanClass;   % Kir
Kir.name = 'Kir';
Kir.ChanFunc = @Chan_Kir_Struyk;
Kir.molecs_needed = [2];
Kir.ion_current_mult = [2 1];
Kir.eta = 1;

Kir.use_GHK = 1;

Kir.fit_params = [];

%Kir.params(1) = 1e-5;       % P_Kir, cm/s
%Kir.params(1) = 0.25*1e-5;       % P_Kir, cm/s
Kir.params(1) = 0.35*1e-5;       % P_Kir, cm/s

%Kir.params(1) = 6e-6;       % P_Kir, cm/s

%Kir.params(1) = 4e-5;       % P_Kir, cm/s
%Kir.params(1) = 2e-6;       % P_Kir, cm/s      combine with Cl = 2e-4 for TC of ~3ms
%Kir.params(1) = 0;       % P_Kir, cm/s

%Kir.params(1) = 1.4e-6;       % P_Kir, cm/s     combine with P_Cl of 2e-4 for TC of ~3ms, and P_Cl of 2e-5 for TC of ~ 12ms

%G_Kir = 3.7;        % ms/cm^2      From Wallinga
%P_Kir = permeability_from_conductance(G_Kir,K_o,K_i,1,Temperature);
%G_Kir = 0.0;     % mS/cm^2/mM    DiFranco using permeability of 1 x 10^-5 cm/s and temp of 298K                    

Kir.params(2) = -19.5;    %delta_V50_S1
Kir.params(3) = 13.7;  %Aslope_S1
Kir.params(4) = 5;  %tauS1_mag
Kir.params(5) = 0;  %tauS1_Vmean
Kir.params(6) = 500;    %tauS1_Vwidth


Kir.params(7) = -200; %delta_V50_S2, mV
Kir.params(8) = 5; %Aslope_S2, mV
Kir.params(9) = 5; %tauS2_mag, ms
Kir.params(10) = 0; %tauS2_Vmean, mV
Kir.params(11) = 500;  %tauS2_Vwidth, mV

% params below are not meant to be fit
Kir.params(12) = 1; %valence   % valence of ion
Kir.params(13) = 1; %S1_act    % -1 for activating, 1 for inactivating
Kir.params(14) = 1; %S1_exp    % exponent when determining final open fraction (such as the 3 in m^3)
Kir.params(15) = -1; %S2_act
Kir.params(16) = 1; %S2_exp

Kir.param_bounds = [1 3 0;
    
                2 20 0;   
                1 3 0;    
                1 3 0;    
                2 20 0;
                1 3 0;
                
                2 20 0;   
                1 3 0;    
                1 3 0;    
                2 20 0;
                1 3 0;
                
                1 1 0;   
                1 1 0;    
                1 1 0;    
                1 1 0;
                1 1 0];    

[Kir.statenames, Kir.vars] = Kir.ChanFunc([],[],[],1,[],[],[],[],[]); 
Kir.vars{3} = [0.496, 0.9999]; 
Chan_base(chan_number) = Kir;
end

            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cl params
if 1
chan_number = chan_number + 1;
Clc1 = ChanClass;   % Clc1
Clc1.name = 'ClC-1';
Clc1.ChanFunc = @Chan_Cl_DiFranco;
Clc1.molecs_needed = [3];
Clc1.ion_current_mult = [3 1];
Clc1.eta = 1;

Clc1.use_GHK = 1;

Clc1.fit_params = [];

%G_Cl_max = 30.0;   % mS/cm^2
%P_Cl_max = permeability_from_conductance(G_Cl_max,Cl_o,Cl_i,1,Temperature);    % cm/s

%Clc1.params(1) = 6e-5;   %Cl permeability, cm/s
Clc1.params(1) = 3e-4;   %Cl permeability, cm/s    combine with Kir = 2e-6 or 1e-6 for TC of ~3ms
%Clc1.params(1) = 4e-4;   %Cl permeability, cm/s    combine with Kir = 2e-6 or 1e-6 for TC of ~3ms
%Clc1.params(1) = 2e-5;   %Cl permeability, cm/s
%Clc1.params(1) = 2e-6;   %Cl permeability, cm/s

Clc1.params(2) = -50;  % V_alpha_Cl_bar, mV
Clc1.params(3) = 41; %  K_alpha_Cl, mV
Clc1.params(4) = 0.03; % alpha_Cl_bar, ms^-1
Clc1.params(5) = -90;  %  V_beta_Cl_bar, mV
Clc1.params(6) = 25;    % K_beta_Cl, mV
Clc1.params(7) = 0.16;       % beta_Cl_bar, 1/ms
Clc1.params(8) = 1.5;     % Cl_m, mM, Effective extracellular Cl concentration due to surface "obstruction".  Only in use if use_GHK = 1.

Clc1.param_bounds = [1 3 0;
    
                2 30 0;   
                2 30 0;
                1 4 0;
                
                2 30 0;    
                2 30 0;    
                1 4 0;
                
                1 4 0];   

[Clc1.statenames, Clc1.vars] = Clc1.ChanFunc([],[],[],1,[],[],[],[],[]); 
Clc1.vars{3} = [0.1027]; 
Chan_base(chan_number) = Clc1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NaK params

% expect net NaK current at rest to be around 300 nA/cm^2 (based on calculation from Clausen 2013 J Gen Physiol, pg 339 value of 0.287 umole/g wet/min for Na current
% and intra volume = 75% of total and cell is cylinder)
if 1
chan_number = chan_number + 1;
NaK = ChanClass;   % NaK ATPase
NaK.name = 'NaK';
NaK.ChanFunc = @Chan_NaK_Wallinga;
NaK.molecs_needed = [1 2];
NaK.ion_current_mult = [1 3;2 -2];
NaK.eta = 0.1;
NaK.conductance_constant = 96485.33;

NaK.use_GHK = 0;    % generally not applicable

NaK.fit_params = [];

%NaK.params(1) = 0; % umol/(cm^2 s), J_NaK
%NaK.params(1) = 207e-6; % umol/(cm^2 s), J_NaK
NaK.params(1) = 40e-6; % umol/(cm^2 s), J_NaK
%NaK.params(1) = 20e-6; % umol/(cm^2 s), J_NaK
NaK.params(2) = 1;  % mM, Km_K;
NaK.params(3) = 13; % mM, Km_Na;
NaK.params(4) = 0.12; % dimensionless, c1 
NaK.params(5) = 0.1;   % dimensionless, c2
NaK.params(6) = 0.04;   % dimensionless, c3
NaK.params(7) = 7;      % dimensionless, c4
NaK.params(8) = 67.3;   % mM, c5


NaK.param_bounds = [1 3 0;
    
                1 3 0;   
                1 3 0;
                
                1 2 0;
                1 3 0;    
                1 3 0;    
                1 3 0;
                1 3 0];   

[NaK.statenames, NaK.vars] = NaK.ChanFunc([],[],[],1,[],[],[],[],[],[]); 
Chan_base(chan_number) = NaK;

%new_PNaK = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NaPIC params
if 1
chan_number = chan_number + 1;
NaPIC = ChanClass;   % NaPIC
NaPIC.name = 'NaPIC';
NaPIC.ChanFunc = @Chan_generic_2HHstates;
NaPIC.molecs_needed = [1];
NaPIC.ion_current_mult = [1 1];
NaPIC.eta = 1;

NaPIC.use_GHK = 1;

NaPIC.fit_params = [];

%G_NaPIC_max = 0.05;   % mS/cm^2
%G_NaPIC_max = 0.025;   % mS/cm^2
%G_NaPIC_max = 0.065;   % mS/cm^2
%G_NaPIC_max = 0.1;   % mS/cm^2
%G_NaPIC_max = 0;   % mS/cm^2
%Na_o = Molecs_base(1).init_extrac;
%Na_i = Molecs_base(1).init_intrac;
%P_NaPIC_max = permeability_from_conductance(G_NaPIC_max,Na_o,Na_i,1,Temperature);    % cm/s
%NaPIC.params(1) = P_NaPIC_max;   %NaPIC permeability, cm/s

%NaPIC.params(1) = 1.8e-8;   %NaPIC permeability, cm/s
NaPIC.params(1) = 0.9e-8;   %NaPIC permeability, cm/s

NaPIC.params(2) = -67.5;    %V50_S1, mV 
NaPIC.params(3) = 3;    %Aslope_S1, mV
NaPIC.params(4) = 5;    %tauS1_mag, ms
NaPIC.params(5) = -50;    %tauS1_Vmean, mV
NaPIC.params(6) = 100;    %tauS1_Vwidth, mV

NaPIC.params(7) = -20;    %V50_S2, mV
NaPIC.params(8) = 3;    %Aslope_S2, mV
NaPIC.params(9) = 5;    %tauS2_mag, ms
NaPIC.params(10) = -50;   %tauS2_Vmean, mV
NaPIC.params(11) = 100;   %tauS2_Vwidth, mV

% params below are not meant to be fit
NaPIC.params(12) = 1;   %valence =    % valence of ion
NaPIC.params(13) = -1;   %S1_act =     % -1 for activating, 1 for inactivating
NaPIC.params(14) = 1;   %S1_exp =     % exponent when determining final open fraction (such as the 3 in m^3)
NaPIC.params(15) = 1;   %S2_act = 
NaPIC.params(16) = 1;   %S2_exp = 

%new_NaPIC_1 = -43;
%new_NaPIC_2 = -53;


NaPIC.param_bounds = [1 3 0;    %1
    
                2 20 0;   %2
                1 3 0;    %3
                1 4 0;    %4
                2 20 0;    %5
                1 8 0;    %6
                
                2 20 0;   %7
                1 3 0;    %8
                1 4 0;    %9
                2 20 0;    %10    
                1 8 0;    %11
                
                1 1 0;   %12
                1 1 0;    %13
                1 1 0;    %14
                1 1 0;    %15
                1 1 0];    %16

            
[NaPIC.statenames, NaPIC.vars] = NaPIC.ChanFunc([],[],[],1,[],[],[],[],[]); 
NaPIC.vars{3} = [0, 1];
Chan_base(chan_number) = NaPIC;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Na of NaKCl
% to account for transporters that don't vary with membrane voltage
if 1
chan_number = chan_number + 1;
NaofNaKCl = ChanClass;   % NaConst
NaofNaKCl.name = 'NaofNaKCl';
NaofNaKCl.ChanFunc = @Chan_NaKCl_cotransport_sat;
NaofNaKCl.molecs_needed = [1 2 3];
NaofNaKCl.ion_current_mult = [1 1];
NaofNaKCl.eta = 1;
NaofNaKCl.conductance_constant = 1;

NaofNaKCl.use_GHK = 0;    % generally not applicable

NaofNaKCl.fit_params = [];

NaofNaKCl.params(1) = 3.7e-8; 
%NaofNaKCl.params(2) = 1e6;  % max effective concentration difference (typical actual conc diff is around 10e6)
NaofNaKCl.params(2) = 1.2e6;  % max effective concentration difference (typical actual conc diff is around 10e6)
NaofNaKCl.params(3) = 7e6;  % level of actual conc difference at which effective conc diff is half of max effective conc diff

NaofNaKCl.param_bounds = [1 3 0];   

[NaofNaKCl.statenames, NaofNaKCl.vars] = NaofNaKCl.ChanFunc([],[],[],1,[],[],[],[],[]); 
Chan_base(chan_number) = NaofNaKCl;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% K of NaKCl
% to account for transporters that don't vary with membrane voltage
if 1
chan_number = chan_number + 1;
KofNaKCl = ChanClass;   % KConst
KofNaKCl.name = 'KofNaKCl';
KofNaKCl.ChanFunc = @Chan_NaKCl_cotransport_sat;
KofNaKCl.molecs_needed = [1 2 3];
KofNaKCl.ion_current_mult = [2 1];
KofNaKCl.eta = 1;
KofNaKCl.conductance_constant = 1;

KofNaKCl.use_GHK = 0;    % generally not applicable

KofNaKCl.fit_params = [];

KofNaKCl.params(1) = NaofNaKCl.params(1); 
KofNaKCl.params(2) = NaofNaKCl.params(2); 
KofNaKCl.params(3) = NaofNaKCl.params(3); 
%KofNaKCl.params(1) = 3.7e-8; 
%KofNaKCl.params(2) = 1e6; 
%KofNaKCl.params(3) = 7e6; 

KofNaKCl.param_bounds = [1 3 0];   

[KofNaKCl.statenames, KofNaKCl.vars] = KofNaKCl.ChanFunc([],[],[],1,[],[],[],[],[]); 
Chan_base(chan_number) = KofNaKCl;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cl of NaKCl
% to account for transporters that don't vary with membrane voltage
if 1
chan_number = chan_number + 1;
ClofNaKCl = ChanClass;   % ClConst
ClofNaKCl.name = 'ClofNaKCl';
ClofNaKCl.ChanFunc = @Chan_NaKCl_cotransport_sat;
ClofNaKCl.molecs_needed = [1 2 3];
ClofNaKCl.ion_current_mult = [3 1];
ClofNaKCl.eta = 1;
ClofNaKCl.conductance_constant = 2;

ClofNaKCl.use_GHK = 0;    % generally not applicable

ClofNaKCl.fit_params = [];

ClofNaKCl.params(1) = NaofNaKCl.params(1); 
ClofNaKCl.params(2) = NaofNaKCl.params(2); 
ClofNaKCl.params(3) = NaofNaKCl.params(3); 
%ClofNaKCl.params(1) = 3.7e-8; 
%ClofNaKCl.params(2) = 1e6; 
%ClofNaKCl.params(3) = 7e6; 

ClofNaKCl.param_bounds = [1 3 0];   

[ClofNaKCl.statenames, ClofNaKCl.vars] = ClofNaKCl.ChanFunc([],[],[],1,[],[],[],[],[]); 
Chan_base(chan_number) = ClofNaKCl;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% K of KCl
% to account for transporters that don't vary with membrane voltage
if 0
chan_number = chan_number + 1;
KofKCl = ChanClass;   % KConst
KofKCl.name = 'KofKCl';
KofKCl.ChanFunc = @Chan_KCl_cotransport;
KofKCl.molecs_needed = [2 3];
KofKCl.ion_current_mult = [2 1];
KofKCl.eta = 1;
KofKCl.conductance_constant = 1;

KofKCl.use_GHK = 0;    % generally not applicable

KofKCl.fit_params = [];

KofKCl.params(1) = 0.05; % uA/cm^2  negative current means Na ions moving from extracellular to intracellular
%NaConst.params(1) = -0.360; % uA/cm^2  negative current means Na ions moving from extracellular to intracellular
%%NaConst.params(1) = -0.0180; % uA/cm^2  negative current means Na ions moving from extracellular to intracellular
%%NaConst.params(1) = -0.0540; % uA/cm^2  negative current means Na ions moving from extracellular to intracellular
%%NaConst.params(1) = -0.2; % uA/cm^2  negative current means Na ions moving from extracellular to intracellular

KofKCl.param_bounds = [1 3 0];   

[KofKCl.statenames, KofKCl.vars] = KofKCl.ChanFunc([],[],[],1,[],[],[],[],[]); 
Chan_base(chan_number) = KofKCl;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cl of KCl
% to account for transporters that don't vary with membrane voltage
if 0
chan_number = chan_number + 1;
ClofKCl = ChanClass;   % ClConst
ClofKCl.name = 'ClofKCl';
ClofKCl.ChanFunc = @Chan_KCl_cotransport;
ClofKCl.molecs_needed = [2 3];
ClofKCl.ion_current_mult = [3 1];
ClofKCl.eta = 1;
ClofKCl.conductance_constant = 1;

ClofKCl.use_GHK = 0;    % generally not applicable

ClofKCl.fit_params = [];

ClofKCl.params(1) = 0.05; % uA/cm^2  negative current means Na ions moving from extracellular to intracellular
%NaConst.params(1) = -0.360; % uA/cm^2  negative current means Na ions moving from extracellular to intracellular
%%NaConst.params(1) = -0.0180; % uA/cm^2  negative current means Na ions moving from extracellular to intracellular
%%NaConst.params(1) = -0.0540; % uA/cm^2  negative current means Na ions moving from extracellular to intracellular
%%NaConst.params(1) = -0.2; % uA/cm^2  negative current means Na ions moving from extracellular to intracellular

ClofKCl.param_bounds = [1 3 0];   

[ClofKCl.statenames, ClofKCl.vars] = ClofKCl.ChanFunc([],[],[],1,[],[],[],[],[]); 
Chan_base(chan_number) = ClofKCl;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Na YASC (yet another sodium current)
% initially intended for gating pore current

if 1
chan_number = chan_number + 1;
NaGatingPore = ChanClass;   % NaGatingPore
NaGatingPore.name = 'GatingPore';
NaGatingPore.ChanFunc = @Chan_generic_2HHstates;
NaGatingPore.molecs_needed = [1];
NaGatingPore.ion_current_mult = [1 1];
NaGatingPore.eta = 1;

NaGatingPore.use_GHK = 1;

NaGatingPore.fit_params = [];

%G_NaGatingPore_max = 0.0002;   % mS/cm^2
%G_NaYasc_max = 0.025;   % mS/cm^2
%G_NaGatingPore_max = 0.00;   % mS/cm^2
%Na_o = Molecs_base(1).init_extrac;
%Na_i = Molecs_base(1).init_intrac;
%P_NaGatingPore_max = permeability_from_conductance(G_NaGatingPore_max,Na_o,Na_i,1,Temperature);    % cm/s
%NaYasc_params(1) = 7e-4;   %NaYasc permeability, cm/s
%NaGatingPore.params(1) = P_NaGatingPore_max;   %NaYasc permeability, cm/s

%NaGatingPore.params(1) = 2.444e-9;   %NaYasc permeability, cm/s
%NaGatingPore.params(1) = 2.444e-8;   %NaYasc permeability, cm/s
NaGatingPore.params(1) = 0;   %NaYasc permeability, cm/s

NaGatingPore.params(2) = -240;    %V50_S1, mV 
NaGatingPore.params(3) = 1;    %Aslope_S1, mV
NaGatingPore.params(4) = 5;    %tauS1_mag, ms
NaGatingPore.params(5) = -50;    %tauS1_Vmean, mV
NaGatingPore.params(6) = 200;    %tauS1_Vwidth, mV

NaGatingPore.params(7) = -60;    %V50_S2, mV
NaGatingPore.params(8) = 15;    %Aslope_S2, mV
NaGatingPore.params(9) = 5;    %tauS2_mag, ms
NaGatingPore.params(10) = -50;   %tauS2_Vmean, mV
NaGatingPore.params(11) = 200;   %tauS2_Vwidth, mV

% params below are not meant to be fit
NaGatingPore.params(12) = 1;   %valence =    % valence of ion
NaGatingPore.params(13) = -1;   %S1_act =     % -1 for activating, 1 for inactivating
NaGatingPore.params(14) = 1;   %S1_exp =     % exponent when determining final open fraction (such as the 3 in m^3)
NaGatingPore.params(15) = 1;   %S2_act = 
NaGatingPore.params(16) = 1;   %S2_exp = 

NaGatingPore.param_bounds = [1 3 0;    %1
    
                2 20 0;   %2
                1 3 0;    %3
                1 4 0;    %4
                2 20 0;    %5
                1 8 0;    %6
                
                2 20 0;   %7
                1 3 0;    %8
                1 4 0;    %9
                2 20 0;    %10    
                1 8 0;    %11
                
                1 1 0;   %12
                1 1 0;    %13
                1 1 0;    %14
                1 1 0;    %15
                1 1 0];    %16

            
[NaGatingPore.statenames, NaGatingPore.vars] = NaGatingPore.ChanFunc([],[],[],1,[],[],[],[],[]); 
NaGatingPore.vars{3} = [1, 1];
Chan_base(chan_number) = NaGatingPore;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Na leak current, for resting potential Na permeability
% might be due to regular NaV channel or NaPIC channel, but separating out so don't need 
% to adjust those voltage dependent channels to make sure there is Na leak near resting potential

if 1
chan_number = chan_number + 1;
NaLeak = ChanClass;   % NaLeak
NaLeak.name = 'NaLeak';
NaLeak.ChanFunc = @Chan_cond_only;
NaLeak.molecs_needed = [1];
NaLeak.ion_current_mult = [1 1];
NaLeak.eta = 1;

NaLeak.use_GHK = 1;

NaLeak.fit_params = [];

%G_NaYasc_max = 0.0002;   % mS/cm^2
%G_NaYasc_max = 0.025;   % mS/cm^2

%G_NaLeak_max = 0.00;   % mS/cm^2
%Na_o = Molecs_base(1).init_extrac;
%Na_i = Molecs_base(1).init_intrac;
%P_NaLeak_max = permeability_from_conductance(G_NaLeak_max,Na_o,Na_i,1,Temperature);    % cm/s

%P_NaLeak_max = 2e-8; % cm/s
%P_NaLeak_max = 5e-9; % cm/s
P_NaLeak_max = 2e-9; % cm/s
NaLeak.params(1) = P_NaLeak_max;   %NaLeak permeability, cm/s

NaLeak.param_bounds = [1 3 0];    %1
            
[NaLeak.statenames, NaLeak.vars] = NaLeak.ChanFunc([],[],[],1,[],[],[],[],[]); 
Chan_base(chan_number) = NaLeak;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pgfile = evalc(['type ' mfilename]);    % to save this file as char vector in .mat file for this run

