classdef ChanClass
    properties
        ChanFunc        % function to be used for channel calculations
        eta = 1;        % density of channel in t-tubules relative to sarcolemma
        use_GHK = 0;
        fit_params = [];
        params
        extra_params = [];
        param_bounds
        statenames
        vars
        name
        molecs_needed = [1 2];   % molecule number, according to Molecs_base.  must be in order expected by ChanFunc.  If Nernst
                                % is used, then 1st molecule listed here will be used.  
        conductance_constant = [];   % An empty value is for normal channels.  A non-empty setting is for those channels/currents
                                                    % such as constant currents (where this value = 1) or Na/K ATPase (where this value = F)
                                                    % where the conductance factor from find_max_cond_factor.m doesn't need ion concentration as input.
        ion_current_mult = [1 1];       % connection between molecule and current direction.  A 2-d matrix with 2 columns.  Each row is a different molecule
                                % 1st column is molecule id number.  2nd column is multiplier.  So if current for this channel is 1 nA, and multiplier
                                % is 1, then 1 nA of that ion will be flowing out of cell.  Positive current is from intracellular space to extracellular
                                % space.  For example, for Na-K ATPase, matrix would look like [1 3; 2 -2], assuming molec 1 is Na and molec 2 is K.
    end
end