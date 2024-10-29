% based on 18c_test
% implementing ability to change [K] extracellular at specific times

% initial parameters

% need to run 'format_full_data' first

save_file_prefix = 'Rich_Kidea_19hypoK_norm_to_low_Struyk_v6b';  % for saving output runs


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% membrane parameter file
%pg_10_18_21_Rich_06_AP_repstim22
%pg_10_18_21_Rich_06_AP_repstim22_lowNaP
%pg_10_18_21_Rich_06_AP_repstim21
%pg_10_18_21_Rich_06_AP_repstim21_lowNaP
%pg_10_18_21_Rich_06_AP_repstim21_verylowNaP
%pg_10_18_21_Rich_06_AP_repstim23

%pg_10_18_21_Rich_06_AP_repstim21_slowinact
%pg_Rich_myo_from_Phil_expmtl
%pg_10_18_21_Rich_06_AP_repstim21_lowNaP_slowinact
%pg_10_18_21_Rich_06_AP_repstim21_lowNaT_slowinact

pg_Rich_Kidea_19hypoK_norm_to_low_Struyk_v6b


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% each data file contains FileData.times, FileData.voltage_traces, FileData.current_traces
% times will be 1-d vector.  traces will have rows for each time, and may have multiple columns for "steps"
% for subtracted data sets (TTX, P/4), store the resultant subtracted values as own file

%filename{1} = '10-18-21_Rich_06_mp';
%filename{2} = 'Phil_22_03_19_F2_inact_P4sub';

% list of protocols to be run.  May be 1 for each filename, or may use filename multiple times
% protocols will contain index for filename, type of input (voltage trace, current trace, something else see below), 
% times to use for simulation, times to use for fit, times to use for zero values
% time to use for step lavels, whether zero values change for each step (for example, initial state before inactivation final step is not uniform),
% steps to run, type of objective function for fitting

% types of input:
% 0 for using a voltage trace (from voltage electrode)
% 1 for a current trace (so fit to voltage), currently a vclamp protocol
% 2 for a step to specified current for a specified
%   duration (see 'input_current_value',
%   'input_current_start', and 'input_current_end'
%   below)
% 3 for using the rig model, so voltage clamp signal is
%   the input.  Voltage clamp signal currently is a sigmoidal function, 
%   specified by step_rise and step_midtime toward bottom of this file
% 4 use sigmoid voltage directly, with init and final
%   voltage determined by voltage trace
% 5 use sigmoid voltage directly, with init and final
%   voltage determined by this file.  Only single
%   voltage change allowed.
% 6 for sequence of linear ramps and holds to new voltage values.   Uses 'sim_voltage_input' variable below.
% 7 for constructed current protocols.  Current values will be interpolated at time points as needed by simulation.  See 'construct_current_protocol.m'.
%   Must load 'ic_mat' file before running sim.

% 
protocol{1}.name = '200 ms stim';       % a label, enter anything
%protocol{1}.steps_to_run = [8 10 14 16];    % optional field.  default = 1;
protocol{1}.filenum = 0;    % 0 for pure simulation
protocol{1}.input = 2;
%protocol{1}.zero_times_after = [86 88];     % used for P/4 and TTX subtracted traces, used to force traces to go to zero at end
                                            % code checks for existence of this field, and subtracts mean current over these times
                                            % remove field for other protocols
%protocol{1}.uniform_initial_steady_state = 1;   % field only necessary if multiple steps_to_run
%protocol{1}.step_label_time = 89;       % field only necessary if multiple steps_to_run
%protocol{1}.ic_mat = ic_mat;
%protocol{1}.forced_time_evals = forced_time_evals;
protocol{1}.obj_func = 0;   % placeholder
protocol{1}.time_to_establish_steady_state = 1e-4;

protocol{1}.Vminit_method = 1;
protocol{1}.RM_orig = -74.6;

% repetitive stim section
if 0
% for 0.1 cm length, 5 segments
current_stim_start = 75841;
current_stim_duration = 200; % ms
current_stim_freq = 20; % Hz
current_stim_num = 1;  % number of current stim pulses
%current_stim_density = 9000;   % nA/cm^2

current_stim_density = 65000;   % nA/cm^2
%current_stim_density = 10000;   % nA/cm^2

%current_stim_density = 450000;   % nA/cm^2
%num_current_steps = 10;
%fraction_of_duration_for_steps = 0.2;
num_current_steps = 1;  % use number > 1 to make it so current doesn't jump up to final value.
fraction_of_duration_for_steps = 0;

this_cs_start = current_stim_start;
input_current_density = [];
delta_csstart = 1000/current_stim_freq;

for ncs = 1:current_stim_num
    % make current rise go in 10 steps, covering 1st 20% of current stim duration    
    cst = this_cs_start;
    for csdf = 1:num_current_steps
        frac = csdf/num_current_steps;
        input_current_density = [input_current_density; cst current_stim_density*frac];
        cst = cst + 1/num_current_steps * current_stim_duration * fraction_of_duration_for_steps;
    end
    % make current fall in 10 steps, covering last 20% of current stim duration
    cst = this_cs_start + current_stim_duration * (1 - fraction_of_duration_for_steps);
    for csdf = 1:num_current_steps
        frac = (num_current_steps-csdf)/num_current_steps;
        cst = cst + 1/num_current_steps * current_stim_duration * fraction_of_duration_for_steps;
        input_current_density = [input_current_density; cst current_stim_density*frac];
    end
    this_cs_start = this_cs_start + delta_csstart;
end
input_current_density_start = 0;  % current at beginning steady state, and at beginning of sim (nA/cm^2)
end

if 1
% for 0.1 cm length, 5 segments
input_current_density = [10000 0;
                20000 0];
input_current_density_start = 0;  % current at beginning steady state, and at beginning of sim (nA/cm^2)
end


protocol{1}.input_current_density = input_current_density;
protocol{1}.input_current_density_start = input_current_density_start;  % current at beginning steady state, and at beginning of sim (nA/cm^2)

if exist('current_stim_freq','var')
    zero_times_AP = [current_stim_start - 19 current_stim_start - 9];
    sim_times_AP = [current_stim_start - 9 this_cs_start + 1500];
    fit_times_AP = [32479.0 32483.5];
else

    %steps_to_fit_misc = (1:9);   % indices, not voltages
    zero_times_AP = [49605 49606];     % in ms, [start end]
    %sim_times_AP = [32470 32500];    % in ms, only 1 pair allowed. [start end]
    sim_times_AP = [49606 49630];    % in ms, only 1 pair allowed. [start end]
    %fit_times = [71 71.5;   % in ms,  fit data between 2 times in each row.
    %                88 90];     % All included times must be subset of sim_times_act.
    fit_times_AP = [49606.9 49610];
    %step_label_misc = 27;  % in ms, time to use for voltage labels
end

protocol{1}.sim_times = sim_times_AP;
protocol{1}.fit_times = fit_times_AP;
protocol{1}.zero_times_before = zero_times_AP;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if input = 2, then include following in protocol:

%protocol{n}.input_current_density = [32479.37 60000;
%                32480.385 0];
%protocol{n}.input_current_density_start = 0;  % current at beginning steady state, and at beginning of sim (nA/cm^2)

% for input_current, first column is time for new current, second column is
% value of current density that starts at that time (in nA/cm^2)
% Times used must be within the 'sim_times' values.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if filenum = 0 (pure sim), include the following:

protocol{1}.created_delta_t = 1000; % ms
protocol{1}.final_t = 6e7;   % ms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if input = 6, simmed voltage sequence, include following:

% first column is times (in ms).  Rest of columns contain voltages and are for protocol steps.
% So a single 'experiment' will only have 2 columns.

% First row always has time of 0, and each other column contains initial
% voltage.

% Can't have instantaneous voltage steps.  Must use a different time if
% wish a voltage change.

% Code will assume last voltage given is held
% until end of time range simmed.

% protocol{n}.sim_voltage_input = [0 -85;
%                                 5 -85;
%                                 5.1 -30];
                    
% multiple step example. All columns must have same number of rows, so very 
% different protocols can't be done with this.  Need to create a different
% matrix (and do a separate run) for a different protocol.
% protocol{n}.sim_voltage_input = [0 -85 -85 -85 -85;
%                                  5 -85 -85 -85 -85;
%                                  5.1 -50 -40 -30 -20];
                


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if input == 4 or 5 (sigmoid voltage input), include following:

%protocol{n}.sig_step_rise = 33; % sharpness of transition
%protocol{n}.sig_step_midtime = 72.93; % midpoint of voltage rise (ms)

% following only used if input = 5
%protocol{n}.sig_step_init = 0;    % mV, only for input_method 5, can be a vector for protocol steps (be sure to include steps_to_run if multiple steps)
%protocol{n}.sig_step_final = 31;   % mV, only for input_method 5, can be a vector for protocol steps (same length as above line)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if different protocols have different initial ion concentrations, can include following
% any concentrations not specified will assume default values stated in pg_ file
% index matches that used for ion in pg_ file

% protocol{n}.conc.modified_molecs = [1 3]
% protocol{n}.conc.Molecs_base(1).init_extrac = 150
% protocol{n}.conc.Molecs_base(1).init_intrac = 15
% protocol{n}.conc.Molecs_base(1).init_tt = 150     % not used

%protocol{1}.conc.modified_molecs = [1];
%protocol{1}.conc.Molecs_base(1).init_extrac = 143;

% protocol{n}.conc.osmolarity = 290
% protocol{n}.conc.zX = -1.6

% etc.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if different protocols turn off or reduce (or increase?) permeability of specific channels, can include following
% typically used if drug blocks a channel
% any permeabilities not specified will assume default values stated in pg_ file
% identify channel by .name field in pg_ file
% .modified_perms is fraction of value specified in pg_ file

% protocol{n}.modified_perm_chans = ["Kdr" "NaK"];
% protocol{n}.modified_perms = [0.01 0.0];
% protocol{1}.modified_perm_chans = ["Kdr"];
% protocol{1}.modified_perms = [0.5];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for adding or subtracting extracellular ions at specific time points
% only works with input methods 1, 2, or 7

% values below are for mM (millimolar) changes for ions
% calculations will assume *initial* extracellular volume
% Thus if conc change below is -3 mM, then the total number of millimoles 
% of ions to be removed will be 3*Volume_ex, with Volume_ex the initial 
% extracellular volume in liters.
% Thus this section is used to specify adding or removing a specific total 
% number of ions.

% This change will be transformed into an *ion current* that occurs over the
% specified time frame.
% Assumes sim starts with no change to ion concentrations specified in pg file.

% 1st column is time to begin change (must be within simtimes)
% 2nd column is the *change* in mM required
% 3rd column is the duration over which to effect the change (in ms)

% Warning:  It is possible to make extracellular concentrations go negative
% if you take out more than exists *at the time the change occurs* in a 
% short time period.  Usually this will cause an integration tolerances
% failure

% Should keep these changes to charge neutral.  Not sure what happens
% otherwise.

protocol{1}.molecs_change.modified_molecs = [1 2];
%protocol{1}.molecs_change.molec(1).mat = [4e7 5 1000000];
%protocol{1}.molecs_change.molec(2).mat = [4e7 -5 1000000];
%protocol{1}.molecs_change.molec(1).mat = [3e7 -12 10000];
%protocol{1}.molecs_change.molec(2).mat = [3e7 12 10000];
%protocol{1}.molecs_change.molec(1).mat = [3e7 4 10000];
%protocol{1}.molecs_change.molec(2).mat = [3e7 -4 10000];
protocol{1}.molecs_change.molec(1).mat = [3e7 3 10000];
protocol{1}.molecs_change.molec(2).mat = [3e7 -3 10000];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set ion being passed for electrode current
% default is K (mol num 2 actually) and Cl (mol number 3)

% protocol{n}.electrode_cation = 1
% protocol{n}.electrode_anion = ?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% can include following:

%protocol{n}.tshift = 0.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
do_fits_or_onerun = 2;  % 0 for neither (for just plotting raw data)
                        % 1 for doing fit
                        % 2 for doing one run with parameters
                        %   The 
                        %   program will use init_guess parameters.  
                        %   Changing set of init_guess params in
                        %   param file is safest way.
                        
use_IC_comps = 0;   % use intracellular compartments.  (TT compartments and TT radial current are not turned on and off by this)
                    % Longitudinal IC current is *not* turned off by this.
                    % Longitudinal TT current *is* turned off by this.
                    % Radial IC current *is* turned off by this.
                    
use_prev_fit_params = 0;
time_resolution_increase = 1;    % Use 1 if match with experimental data
                            % higher numbers mean to add N-1 extra
                            % simulation output points.  So 10 means 9
                            % extra time points for each real time point.
                            

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% output control
plotstep = 1;

do_basic_graphs = 1;
do_medium_graphs = 1;
do_advanced_output = 0; % advanced output includes capacitance calc,
                        % I-V plot, current peak calcs, etc.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fit algorithm control
fit_algorithm = 2;   % 1 for lsqcurvefit (uses gradients), 2 for patternsearch
MaxFunctionEvaluations = 320;

% following used only for lsqcurvefit
MaxIterations = 100;
FiniteDifferenceStepSize = 1e-2;
FiniteDifferenceType = 'central';    %choices are 'forward' or 'central'
StepTolerance = 1e-5;
FunctionTolerance = 1e-5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sim algorithm control
%RelTol = 1e-4;
%AbsTol = 1e-7;
RelTol = 1e-8;
AbsTol = 1e-11;
%RelTol = 1e-11;
%AbsTol = 1e-14;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rig stuff
step_rise = 9;
step_midtime = 23.51;

Re1 = 50;    % MOhms
Re2 = 50;    % MOhms, not used atm (because using current source with R0)
%tau_mu = 0.02;  % ms
tau_mu = .02;  % ms
C_in = 0.004;   % nF
%C_n = 0.002;    % nF
C_n = 0.0039;    % nF
mu_gain = 10000;  % gain
R0 = 10;     %MOhms
                            
pfile = evalc(['type ' mfilename]); % to save this file as char vector in .mat file for this run




