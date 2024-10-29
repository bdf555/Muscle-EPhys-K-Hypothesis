function [pred_vals, pred_vals_sameunits, outstruct, state_num] = track_vol_calcs(state_num_init, ode_result, sP, outstruct, nprot)

% rearrange variables so nsegs is 1st dimension, shells is 2nd, so functions work properly
% and store molecule numbers
state_num = state_num_init;
pred_vals = ode_result;
pred_vals_sameunits = ode_result;
for nms = 1:sP.num_molecs
    state_num = state_num + 1;
    outstruct.prot{nprot}.molecs_t_rsh_n{nms} = ode_result{state_num};
    rmolect{nms} = permute(ode_result{state_num},[3,4,1,2]);
    state_num = state_num + 1;
    outstruct.prot{nprot}.molecs_i_n{nms} = ode_result{state_num};
    rmoleci{nms} = permute(ode_result{state_num},[3,1,2]);
    state_num = state_num + 1;  % don't need extracellular concentrations, so skip
    outstruct.prot{nprot}.molecs_e_n{nms} = ode_result{state_num};
    rmolece{nms} = ode_result{state_num};
end

[volume_t_rsh, volume_i, volume_e, ~, ~, ~] = Volume_from_Ion_nums(sP, rmolect, rmoleci, rmolece);    

[Vmcalc, Vtcalc, Cmolect, Cmoleci, Cmolece] = Volt_and_Conc_from_Ion_Nums(sP, rmolect, rmoleci, rmolece, volume_t_rsh, volume_i, volume_e);

pred_vals{1} = permute(Vmcalc,[2,3,1]); 
pred_vals_sameunits{1} = pred_vals{1};
pred_vals{4} = permute(Vtcalc,[3,4,1,2]); 
pred_vals_sameunits{4} = pred_vals{4};

pred_vals{5} = permute(volume_t_rsh,[3,4,1,2]); 
pred_vals_sameunits{5} = pred_vals{5}/1000;
pred_vals{6} = permute(volume_i,[2,3,1]); 
pred_vals_sameunits{6} = pred_vals{6}/1000;
pred_vals{7} = volume_e; 
pred_vals_sameunits{7} = pred_vals{7}/1000;

% rearrange variables so normal dimension order (times, steps, segs, shells)
% and store concentrations
state_num = state_num_init;
for nms = 1:sP.num_molecs
    state_num = state_num + 1;
    if sP.num_shells > 0
        pred_vals{state_num} = permute(Cmolect{nms},[3,4,1,2]);
        outstruct.prot{nprot}.molecs_t_rsh{nms} = pred_vals{state_num};
    else
        outstruct.prot{nprot}.molecs_t_rsh{nms} = [];
    end
    state_num = state_num + 1;
    pred_vals{state_num} = permute(Cmoleci{nms},[2,3,1]);
    outstruct.prot{nprot}.molecs_i{nms} = pred_vals{state_num};
    state_num = state_num + 1;  
    pred_vals{state_num} = Cmolece{nms};
    outstruct.prot{nprot}.molecs_e{nms} = pred_vals{state_num};
end
