
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
cell_length = 1; % cm


variable_intracellular_volume = 1;     % intracellular volume is allowed to change to keep osmolarity constant. Was called 'track_volumes' before.
infinite_extracellular_volume = 0;  % if 1, then renders extrac_volume_multiple meaningless.  This value is not used in ode version (track time course)
                                    % since the code has no provision for infinite extracellular volume.  This value is used in the steady state solver version.
track_volumes = variable_intracellular_volume;

extrac_volume_multiple = 0.5;     % extracellular volume = intracellular volume times this multiple
                                % interstitial muscle volume is about 25% (0.25) of the total muscle volume (Kandarian, 1991, J Appl Physiol
use_electrical_drift = 0;   % when currents flow along tube, they are carried by ions.  So ion concentrations can change this way.
                            % Must use if t-tubules are on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% output focus
plotseg = voltage_probe_segment;
plotshell = voltage_probe_shell;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Temperature = 305;			% K

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
use_RM_orig = 0;        % This is not used for steady state version, since must have initial RM_orig guess.
RM_orig = -60;          % an initial guess for RM (if variable_intracellular_volume) or the starting point of RM
                        % If use_RM_orig = 0, then initial guess for time course version is set to Nernst for K based on initial values

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if track_volumes
    osmolarity = 300;   % mOsm
    zX = -1.6477;    % valence of impermeant intracellular anion
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

vole = voli*extrac_volume_multiple;


if 1    % no ion shifting
Molecs_base(1) = MolecClass;    % Na
Molecs_base(1).init_extrac = 140;    % normal
Molecs_base(1).init_intrac = 12.3;

Molecs_base(1).init_tt = 149;
Molecs_base(1).name = 'Na';
end

if 1    % no ion shifting
Molecs_base(2) = MolecClass;    % K
Molecs_base(2).init_extrac = 10;     % normal
Molecs_base(2).init_intrac = 95;

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
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fittable params


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% general parameters

fit_gen_params = [];

% Parameter 1, uF/cm^2, membrane capacitance per area
                    % area is surface (if no TT) or total (if TT)
gen_params(1)= 0.9;   % specific capacitance, uF/cm^2
                                

% Parameter 2, mS/cm^2, series conductance
gen_params(2) = 100000000;    % series conductance, ms/cm^2
                                
% Parameter 3, mS/cm, intracellular medium conductivity
gen_params(3) = 18;    % intracellular medium conductivity, mS/cm
                            
% Parameter 4, conductivity t-tubule lumen, mS/cm
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

P_KDr_max = 2.864e-04;
Kdr.params(1) = P_KDr_max;   %Kdr permeability, cm/s

Kdr.params(2) = -47.001;  % V_n_bar, mV
Kdr.params(3) = 4.192; %  K_alpha_n, mV
Kdr.params(4) = 17.42;  %  K_beta_n, mV
Kdr.params(5) = 0.013621;    % alpha_n_bar, 1/(ms*mV)
Kdr.params(6) = 0.083997;       % beta_n_bar, 1/ms

% slow inact
Kdr.params(7) = -40;    % V50, mV
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

% Nav.statenames contains names of state variables.  = {'m', 'h', 'S'}
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

Kir.params(1) = 1e-5;       % P_Kir, cm/s

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

Clc1.params(1) = 3e-4;   %Cl permeability, cm/s    combine with Kir = 2e-6 or 1e-6 for TC of ~3ms

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

%NaK.params(1) = 207e-6; % umol/(cm^2 s), J_NaK
NaK.params(1) = 40e-6; % umol/(cm^2 s), J_NaK
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

ClofNaKCl.param_bounds = [1 3 0];   

[ClofNaKCl.statenames, ClofNaKCl.vars] = ClofNaKCl.ChanFunc([],[],[],1,[],[],[],[],[]); 
Chan_base(chan_number) = ClofNaKCl;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gating pore current

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


P_NaLeak_max = 2e-9; % cm/s

NaLeak.params(1) = P_NaLeak_max;   %NaLeak permeability, cm/s

NaLeak.param_bounds = [1 3 0];    %1
            
[NaLeak.statenames, NaLeak.vars] = NaLeak.ChanFunc([],[],[],1,[],[],[],[],[]); 
Chan_base(chan_number) = NaLeak;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pgfile = evalc(['type ' mfilename]);    % to save this file as char vector in .mat file for this run

