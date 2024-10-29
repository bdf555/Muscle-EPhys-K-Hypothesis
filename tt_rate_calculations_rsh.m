% radial cable model based on Wallinga et. al. (1999) Eur. Biophys. J.
% with corrected K diffusion scaling

function [dVt_dt, t_memb_plus_elec_current, tdiffusion_molrate, tmembrane_current, tdiffusion_tt_to_ec_molrate, telectrical_current_flowdown] = tt_rate_calculations_rsh(Vabs, Vt_rsh, It_curr_rsh, sP, Vm, dVmdt, IT_partial, molecs_t_rsh, molecs_e, osmo_t_mobile, osmo_e_mobile)
% Vt_rsh and current inputs are matrices, nsegs rows, nshells columns.  

Itnet = 0;
for nits = 1:sP.num_chans
    Itnet = Itnet + It_curr_rsh{nits};
end

if sP.track_volumes
    dVt_dt = 0;
else
    if ~sP.use_IC_comps


        if 1

        IT = IT_partial + sP.C.*dVmdt./sP.g_S./sP.ra;

        I_lumen_up = (circshift(Vt_rsh,-1,2) - Vt_rsh).*sP.gT';
        I_lumen_up(:,end) = IT;

        I_lumen_down = (circshift(Vt_rsh,1,2) - Vt_rsh) .* circshift(sP.gT',1,2);
        I_lumen_down(:,1) = 0;

        dVt_dt = 1./sP.C_t.*(I_lumen_up + I_lumen_down - Itnet);

        end
    else

        if sP.try_small_matrix
            VT = Vabs - Vt_rsh(:);
            dVt_dt_temp = 1./sP.C_t(:).*(sP.TT*VT - Itnet(:));
            dVt_dt = reshape(dVt_dt_temp,[sP.num_segs, sP.num_shells]);
        else
            dVt_dt = 0;
        end

    end
end
% track molec concentrations in t-tubule

tmembrane_current = [];
tdiffusion_molrate = [];
t_memb_plus_elec_current = [];
telectrical_current_flowdown = [];

if sP.use_electrical_drift
    
    ns = sP.num_segs;
    nsh = sP.num_shells;

    Vabst = Vm - Vt_rsh;    % assumes extracellular voltage is zero
    Vabst = [Vabst zeros(ns,1)];
    diffVt = diff(Vabst,1,2);   % positive when more distal shells are more positive

    diffidx = (diffVt>=0) + ones(ns,nsh).*[1:nsh];
    hlpr = [1:ns]'.*ones(ns,nsh);
    diffidxln = sub2ind([ns,nsh+1],hlpr,diffidx);

    %Tconct = Na_t_rsh + K_t_rsh + Cl_t_rsh;
    Tconce = ones(ns,1).*(osmo_e_mobile);
    Tconct = [osmo_t_mobile Tconce];
    Tconctuse = Tconct(diffidxln);

    T_lumen_cond = [sP.gT(1:end-1)' 1./sP.ra];
end

for nms = 1:sP.num_molecs
    
    % membrane current
    tmembrane_current{nms} = 0;
    
    for nchs = 1:size(sP.molec_chan{nms},1)
        thischan = sP.molec_chan{nms}(nchs,1);
        thismult = sP.molec_chan{nms}(nchs,2);
        tmembrane_current{nms} = tmembrane_current{nms} + It_curr_rsh{thischan}.*thismult;
    end
    
    % diffusion
    diffmolt = diff(molecs_t_rsh{nms},1,2);
    diffmoltup = [diffmolt zeros(sP.num_segs,1)];
    diffmoltdown = [zeros(sP.num_segs,1) -diffmolt];

    if sP.use_old_wrong_diffusion
        tdiffusion_molrate{nms} = sP.molec_diff_mult{nms} .* sP.volume_tt_compartment_rsh .* 10^3 .* (diffmoltup + diffmoltdown);
        tdiffusion_tt_to_ec_molrate{nms} = -sP.molec_diff_mult{nms}(:,end) .* sP.volume_tt_compartment_rsh(:,end) .* 10^3 .* (molecs_e{nms} - molecs_t_rsh{nms}(:,end));
        tdiffusion_molrate{nms}(:,end) = tdiffusion_molrate{nms}(:,end) - tdiffusion_tt_to_ec_molrate{nms};
    else
        tdiffusion_molrate{nms} = sP.molec_diff_mult{nms} .* sP.volume_tt_compartment_rsh .* 10^3 .* (sP.radius_tt_for_diff(:,2:end) .* diffmoltup + sP.radius_tt_for_diff(:,1:end-1) .* diffmoltdown);  % nanomoles/ms
        tdiffusion_tt_to_ec_molrate{nms} = -sP.molec_diff_mult{nms}(:,end) .* sP.volume_tt_compartment_rsh(:,end) .* 10^3 .* sP.radius_tt_for_diff(:,end) .* (molecs_e{nms} - molecs_t_rsh{nms}(:,end));  % nanomoles/ms
        tdiffusion_molrate{nms}(:,end) = tdiffusion_molrate{nms}(:,end) - tdiffusion_tt_to_ec_molrate{nms};
    end
    
    t_memb_plus_elec_current{nms} = tmembrane_current{nms};    % nA
    
    % radial, electrical drift
    if sP.use_electrical_drift
        molec_t_all = [molecs_t_rsh{nms} ones(ns,1).*molecs_e{nms}];
        molecuse = molec_t_all(diffidxln);

        telectrical_current_flowdown{nms} = T_lumen_cond.*molecuse.*diffVt./Tconctuse;
        telectrical_current_flowup = [zeros(ns,1) -telectrical_current_flowdown{nms}(:,1:end-1)];
        telectrical_current = telectrical_current_flowdown{nms} + telectrical_current_flowup;   % nA
        
        t_memb_plus_elec_current{nms} = t_memb_plus_elec_current{nms} + telectrical_current;   % nA
    end
    
end
    
