function [volume_t_rsh, volume_i, volume_e, osmo_t_mobile, osmo_i_mobile, osmo_e_mobile] = Volume_from_Ion_nums(sP, molecs_t_rsh_n, molecs_i_n, molecs_e_n)

mobile_i_n = 0;
mobile_e_n = 0;
for nms = 1:sP.num_molecs
    mobile_i_n = mobile_i_n + molecs_i_n{nms};
    mobile_e_n = mobile_e_n + molecs_e_n{nms};
end
total_i_n = mobile_i_n + sP.fixed_intra_anion_n + sP.neut_osmo_intrac_n;
total_e_n = mobile_e_n + sP.imperm_extra_anion_n + sP.neut_osmo_extrac_n;

volume_i = total_i_n ./ sP.osmo;    % uL
volume_e = total_e_n ./ sP.osmo;

osmo_i_mobile = mobile_i_n ./ volume_i; % mM
osmo_e_mobile = mobile_e_n ./ volume_e; % mM

if sP.num_shells > 0
    mobile_t_n = 0;
    for nms = 1:sP.num_molecs
        mobile_t_n = mobile_t_n + molecs_t_rsh_n{nms};
    end
    total_t_n = mobile_t_n + sP.imperm_t_anion_n + sP.neut_osmo_t_n;
    volume_t_rsh = total_t_n ./ sP.osmo;
    
    osmo_t_mobile = mobile_t_n ./ volume_t_rsh;
else
    volume_t_rsh = [];
    osmo_t_mobile = [];
end


