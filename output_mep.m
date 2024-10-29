

numprot = fitData.num_protocols;
num_molecs = length(Molecs_base);
num_chans = length(uIpms.Chan_base);

for nprot = 1:numprot
    plot_steps{nprot} = fitData.prot{nprot}.steps_to_run;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    if do_basic_graphs
    
        if fitData.prot{nprot}.filenum ~= 0
    
            figure(1)
            subplot(1,numprot,nprot)
            plot(FileData{fitData.prot{nprot}.filenum}.times, FileData{fitData.prot{nprot}.filenum}.voltage_traces(:,plot_steps{nprot}))
            title('all full voltage traces')

            if ~isempty(FileData{fitData.prot{nprot}.filenum}.current_traces)
                figure(2)
                subplot(1,numprot,nprot)
                plot(FileData{fitData.prot{nprot}.filenum}.times, FileData{fitData.prot{nprot}.filenum}.current_traces(:,plot_steps{nprot}))
                title('all full current traces')
            end
            
        end


        if do_fits_or_onerun ~= 0
            figure(8)
            subplot(1,numprot,nprot)
            vmplot = squeeze(outstruct.prot{nprot}.Vm_pred(:,plotstep,:));
            plot(fitData.prot{nprot}.times_to_sim, vmplot)

            if uIpms.add_capac_coupling
                hold on
                vmplotcc = squeeze(outstruct.prot{nprot}.Vmeas(:,plotstep,:));
                set(gca,'ColorOrderIndex',1)
                plot(fitData.prot{nprot}.times_to_sim, vmplotcc, ':')
                title(['Vm (surface membrane) (solid) and Vmeas (capac coupl) (dotted), each segment, step ' num2str(plotstep)])
                hold off
            else
                title(['Vm (surface membrane), each segment, step ' num2str(plotstep)])
            end
                

            figure(9)
            if nprot == 1
                clf
            end
            subplot(1,numprot,nprot)
            hold on
            num_lines = length(fitData.prot{nprot}.steps_to_run);
            cc = lines(num_lines);
            for i = 1:num_lines
                if fitData.prot{nprot}.filenum ~= 0
                    plot(fitData.prot{nprot}.times_to_sim, fitData.prot{nprot}.volt_traces_to_sim(:,i), ':','Color',cc(i,:));
                end
                plot(fitData.prot{nprot}.times_to_sim/60000, outstruct.prot{nprot}.Vm_pred(:,i,fitData.uIpms.gen.voltage_probe_segment), 'Color',cc(i,:))
                if uIpms.add_capac_coupling
                    plot(fitData.prot{nprot}.times_to_sim, outstruct.prot{nprot}.Vmeas(:,i,fitData.uIpms.gen.voltage_probe_segment), '--', 'Color',cc(i,:))
                end

            end
            hold off
            if ~uIpms.add_capac_coupling
                if fitData.prot{nprot}.filenum ~= 0
                    title(['Vm (surface membrane), expmt vs pred, each step, voltage probe seg'])
                else
                    title(['Vm (surface membrane), pred, each step, voltage probe seg'])
                end
            else
                if fitData.prot{nprot}.filenum ~= 0
                    title(['Vm (surface membrane), expmt(dotted), pred(solid), withcapcoupl(dashed), each step, voltage probe seg'])
                else
                    title(['Vm (surface membrane), pred(solid), withcapcoupl(dashed), each step, voltage probe seg'])
                end
            end


            figure(10)
            subplot(1,numprot,nprot)
            implot = squeeze(outstruct.prot{nprot}.ImTot(:,plotstep,:));
            plot(fitData.prot{nprot}.times_to_sim, implot, fitData.prot{nprot}.times_to_sim, outstruct.prot{nprot}.Itot_m_pred(:,plotstep),'.')
            title(['Im total (surface membrane plus Ra), each segment and total, step ' num2str(plotstep)])

            figure(11)
            if nprot == 1
                clf
            end
            subplot(1,numprot,nprot)
            hold on
            num_lines = length(fitData.prot{nprot}.steps_to_run);
            cc = lines(num_lines);
            for i = 1:num_lines
                if fitData.prot{nprot}.filenum ~= 0
                    plot(fitData.prot{nprot}.times_to_sim, fitData.prot{nprot}.currents_to_sim(:,i), '--','Color',cc(i,:));
                end
                if ~isfield(fitData.prot{nprot}, 'zero_times_after')
                    plot(fitData.prot{nprot}.times_to_sim, outstruct.prot{nprot}.Itot_m_pred(:,i),'Color',cc(i,:))
                    if fitData.prot{nprot}.filenum ~= 0
                        titlestr = 'I memb total (surface membrane plus Ra) experimental vs. predicted, each step';
                    else
                        titlestr = 'I memb total (surface membrane plus Ra) predicted, each step';
                    end
                else
                    
                    Ifittot_pred = 0;
                    for chanids_to_use = fitData.uIpms.chan_ids_for_plot_and_fit
                        name_num = 1;
                        if fitData.uIpms.fit.num_shells > 0
                            Ifittot_pred = Ifittot_pred + outstruct.prot{nprot}.Itot_m_chan{chanids_to_use}(:,i) + outstruct.prot{nprot}.Itot_t_chan{chanids_to_use}(:,i);
                        else
                            Ifittot_pred = Ifittot_pred + outstruct.prot{nprot}.Itot_m_chan{chanids_to_use}(:,i);
                        end
                        if name_num == 1
                            chan_names = Chan_base(chanids_to_use).name;
                        else
                            chan_names = [chan_names ' + ' Chan_base(chanids_to_use).name];
                        end
                    end
                    
                    if fitData.uIpms.fit.num_shells > 0
                        plot(fitData.prot{nprot}.times_to_sim, Ifittot_pred,'Color',cc(i,:))
                    else
                        plot(fitData.prot{nprot}.times_to_sim, Ifittot_pred,'Color',cc(i,:))
                    end
                    if fitData.prot{nprot}.filenum ~= 0
                        titlestr = [chan_names ' current, memb total (surface membrane plus Ra) experimental vs. predicted, each step'];
                    else
                        titlestr = [chan_names ' current, memb total (surface membrane plus Ra) predicted, each step'];
                    end
                end
            end
            hold off
            title(titlestr)
            
            gS_plot = all_params(2).*fitData.G_to_g;
            vs = (outstruct.prot{nprot}.ImTot - outstruct.prot{nprot}.ImT)./gS_plot;

            figure(12)
            if nprot == 1
                clf
            end
            subplot(1,numprot,nprot)
            hold on
            num_lines = length(fitData.prot{nprot}.steps_to_run);
            cc = lines(num_lines);
            for i = 1:num_lines
                if fitData.prot{nprot}.filenum ~= 0
                    plot(fitData.prot{nprot}.times_to_sim, fitData.prot{nprot}.volt_traces_to_sim(:,i), ':','Color',cc(i,:));
                end
                
                V_pred = vs(:,i,fitData.uIpms.gen.voltage_probe_segment) + outstruct.prot{nprot}.Vm_pred(:,i,fitData.uIpms.gen.voltage_probe_segment);
                plot(fitData.prot{nprot}.times_to_sim, V_pred, 'Color',cc(i,:))
            end
            hold off
            if fitData.prot{nprot}.filenum ~= 0
                title(['Vabs of IC surface compart, expmt vs. predicted, voltage probe seg, each step'])
            else
                title(['Vabs of IC surface compart, predicted, voltage probe seg, each step'])
            end

            
            figure(13)
            if nprot == 1
                clf
            end
            subplot(1,numprot,nprot)
            hold on
            for i = 1:num_lines
                if fitData.prot{nprot}.filenum ~= 0
                    plot(fitData.prot{nprot}.times_to_sim, fitData.prot{nprot}.volt_traces_to_sim(:,i), ':','Color',cc(i,:));
                end
                plot(fitData.prot{nprot}.times_to_sim, vs(:,i,fitData.uIpms.gen.voltage_probe_segment), 'Color',cc(i,:))
            end
            hold off
            if fitData.prot{nprot}.filenum ~= 0
                title(['Vacross series resist, expmt Vabs vs. pred across Rs, voltage probe seg, each step'])
            else
                title(['Vacross series resist, pred across Rs, voltage probe seg, each step'])
            end

            
            figure(14)
            subplot(1,numprot,nprot)
            v_int_vs_ground = outstruct.prot{nprot}.Vm_pred + vs;
            vplot = squeeze(v_int_vs_ground(:,plotstep,plotseg));
            if fitData.prot{nprot}.filenum ~= 0
                plot(fitData.prot{nprot}.times_to_sim, fitData.prot{nprot}.volt_traces_to_sim(:,plotstep), ':', fitData.prot{nprot}.times_to_sim, vplot);
                title(['Vabs of IC surface compart, expmt and predicted, step ' num2str(plotstep) ', segment ' num2str(plotseg)])
            else
                plot(fitData.prot{nprot}.times_to_sim, vplot);
                title(['Vabs of IC surface compart, predicted, step ' num2str(plotstep) ', segment ' num2str(plotseg)])
            end

            if fitData.uIpms.fit.use_IC_comps
                figure(15)
                subplot(1,numprot,nprot)
                if fitData.prot{nprot}.filenum ~= 0
                    plot(fitData.prot{nprot}.times_to_sim, fitData.prot{nprot}.volt_traces_to_sim(:,plotstep), ':', fitData.prot{nprot}.times_to_sim, outstruct.prot{nprot}.Vabs(:, plotstep, plotseg, plotshell));
                    title(['Vabs of IC compart, expmt and pred, step ' num2str(plotstep) ', seg ' num2str(plotseg) ', shell ' num2str(plotshell)])
                else
                    plot(fitData.prot{nprot}.times_to_sim, outstruct.prot{nprot}.Vabs(:, plotstep, plotseg, plotshell));
                    title(['Vabs of IC compart, pred, step ' num2str(plotstep) ', seg ' num2str(plotseg) ', shell ' num2str(plotshell)])
                end
            end
            

            if fitData.uIpms.fit.num_shells > 0

                figure(21)
                subplot(1,numprot,nprot)
                vtplot = squeeze(outstruct.prot{nprot}.Vt_pred(:,plotstep,plotseg,:));
                plot(fitData.prot{nprot}.times_to_sim, vtplot)
                title(['Voltage across TT memb, all shells, step ' num2str(plotstep) ', segment ' num2str(plotseg)])

                figure(24)
                subplot(1,numprot,nprot)

                vt_vs_ground = zeros(size(outstruct.prot{nprot}.Vt_pred));
                for j = 1:fitData.uIpms.gen.num_segments
                    vt_vs_ground(:,:,j,:) = v_int_vs_ground(:,:,j) - outstruct.prot{nprot}.Vt_pred(:,:,j,:);
                end
                vtvg = squeeze(vt_vs_ground(:,plotstep,plotseg,:));
                plot(fitData.prot{nprot}.times_to_sim, vtvg)
                title(['Vabs of TT compart, predicted, all shells, step ' num2str(plotstep) ', segment ' num2str(plotseg)])

            end
        end
        
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if do_medium_graphs && do_fits_or_onerun ~= 0

        fignum = 100;
        figure(fignum)
        subplot(1,numprot,nprot)
        plot(fitData.prot{nprot}.times_to_sim, outstruct.prot{nprot}.Itot_m_C(:,plotstep))
        title(['Total cap current across surface membrane, step ' num2str(plotstep)])
        
        if fitData.uIpms.fit.num_shells > 0
            fignum = fignum + 1;
            figure(fignum);
            subplot(1,numprot,nprot)
            plot(fitData.prot{nprot}.times_to_sim, outstruct.prot{nprot}.Itot_m_T(:,plotstep))
            title(['Current due to voltage across Ra (not due to diffusion), step ' num2str(plotstep)])
        end
        
        for ncs = 1:num_chans
            fignum = fignum + 1;
            figure(fignum);
            subplot(1,numprot,nprot)
            plot(fitData.prot{nprot}.times_to_sim, outstruct.prot{nprot}.Itot_m_chan{ncs}(:,plotstep))
            title(['Total ' Chan_base(ncs).name ' current across surface membrane, step ' num2str(plotstep)])
        end

        show_vps_currents = 1;
        
        if show_vps_currents
            
            showseg = fitData.uIpms.gen.voltage_probe_segment;

            net_vps_curr = 0;
            
            t_of_interest = 16959.38;
            tvps_idx = find(fitData.prot{nprot}.times_to_sim > t_of_interest,1);
            
            fignum = fignum + 1;
            figure(fignum);
            subplot(1,numprot,nprot)
            plot(fitData.prot{nprot}.times_to_sim, outstruct.prot{nprot}.ImC(:,plotstep,showseg))
            title(['Cap current across surface membrane, step ' num2str(plotstep)  ', segment ' num2str(showseg)])
            disp(['Surf memb cap current at time of interest = ' num2str( outstruct.prot{nprot}.ImC(tvps_idx,plotstep,showseg) ) ]);
        
            if fitData.uIpms.fit.num_shells > 0

                fignum = fignum + 1;
                figure(fignum);
                subplot(1,numprot,nprot)
                ttCcurr = squeeze(outstruct.prot{nprot}.ItC(:,plotstep,showseg,:));
                ttCcurr = sum(ttCcurr,2);
                plot(fitData.prot{nprot}.times_to_sim, ttCcurr)
                title(['Total cap current across tt membrane, step ' num2str(plotstep) ', segment ' num2str(showseg)])
                disp(['TT memb cap current at time of interest = ' num2str( ttCcurr(tvps_idx) ) ]);

            end
        
            for ncs = 1:num_chans
                fignum = fignum + 1;
                figure(fignum);
                subplot(1,numprot,nprot)
                plot(fitData.prot{nprot}.times_to_sim, outstruct.prot{nprot}.ImChan{ncs}(:,plotstep,showseg))
                title([Chan_base(ncs).name ' current across surface membrane, step ' num2str(plotstep) ', segment ' num2str(showseg)])
                disp([Chan_base(ncs).name ' current across surf memb at time of interest = ' num2str( outstruct.prot{nprot}.ImChan{ncs}(tvps_idx,plotstep,showseg) ) ]);
                net_vps_curr = net_vps_curr + outstruct.prot{nprot}.ImChan{ncs}(:,plotstep,showseg);
                
                if fitData.uIpms.fit.num_shells > 0
                    fignum = fignum + 1;
                    figure(fignum);
                    subplot(1,numprot,nprot)
                    ttChancurr = squeeze(outstruct.prot{nprot}.ItChan{ncs}(:,plotstep,showseg,:));
                    ttChancurr = sum(ttChancurr,2);
                    plot(fitData.prot{nprot}.times_to_sim, ttChancurr)
                    title(['Total ' Chan_base(ncs).name ' current across tt membrane, step ' num2str(plotstep) ', segment ' num2str(showseg)])
                    disp([Chan_base(ncs).name 'current across tt membrane at time of interest = ' num2str( ttChancurr(tvps_idx) ) ]);
                    net_vps_curr = net_vps_curr + ttChancurr;
                end

            end
            
            if fitData.uIpms.gen.num_segments > 1
                seg_conductance = all_params(3) .* 1000 .* pi .* fitData.uIpms.gen.cell_radius.^2 ./ fitData.deltax;    % uS
                long_curr_left = 0;
                long_curr_right = 0;
                if  showseg > 1
                    long_curr_left = (outstruct.prot{nprot}.Vm_pred(:,i,showseg) - outstruct.prot{nprot}.Vm_pred(:,i,showseg-1) ) .* seg_conductance; %nA
                end
                if showseg < fitData.uIpms.gen.num_segments
                    long_curr_right = (outstruct.prot{nprot}.Vm_pred(:,i,showseg) - outstruct.prot{nprot}.Vm_pred(:,i,showseg+1) ) .* seg_conductance; %nA
                end
                net_long_curr = long_curr_left + long_curr_right;
                
                fignum = fignum + 1;
                figure(fignum);
                subplot(1,numprot,nprot)
                plot(fitData.prot{nprot}.times_to_sim, net_long_curr)
                title(['Net longitudinal current out of segment ' num2str(showseg)])
                disp(['Net long current out of segment ' num2str(showseg) ' at time of interest = ' num2str(net_long_curr(tvps_idx))]);
                net_vps_curr = net_vps_curr + net_long_curr;
                
            end

            fignum = fignum + 1;
            figure(fignum);
            subplot(1,numprot,nprot)
            plot(fitData.prot{nprot}.times_to_sim, net_vps_curr)
            title(['Net current (non-capac) out of segment ' num2str(showseg)])
            disp(['Net current (non-capac) out of segment ' num2str(showseg) ' at time of interest = ' num2str(net_vps_curr(tvps_idx))]);
            
            fignum = fignum + 1;
            figure(fignum);
            subplot(1,numprot,nprot)
            plot(fitData.prot{nprot}.times_to_sim, outstruct.prot{nprot}.Vm_pred(:,plotstep,showseg));
            title(['Vm (surface membrane), step ' num2str(plotstep) ', segment ' num2str(showseg)])
            
            
        end
        
        
        fignum = 199;
        for nms = 1:num_molecs
            
            Molrate_current{nms} = 0;

            for nchs = 1:size(uIpms.molec_chan{nms},1)
                thischan = uIpms.molec_chan{nms}(nchs,1);
                thismult = uIpms.molec_chan{nms}(nchs,2);
                Molrate_current{nms} = Molrate_current{nms} + thismult.*outstruct.prot{nprot}.Itot_m_chan{thischan}(:,plotstep);
                if fitData.uIpms.fit.num_shells > 0
                    Molrate_current{nms} = Molrate_current{nms} + thismult.*outstruct.prot{nprot}.Itot_t_chan{thischan}(:,plotstep);
                end
            end
            
            if length(fitData.prot{nprot}.times_to_sim) == length(Molrate_current{nms})
                fignum = fignum + 1;
                figure(fignum)
                subplot(1,numprot,nprot)

                plot(fitData.prot{nprot}.times_to_sim, Molrate_current{nms})
                title(['Net ' Molecs_base(nms).name ' current (all channels) across surface membrane, step ' num2str(plotstep)])
            end
        end
            
        fignum = 299;
        for ncs = 1:num_chans
            num_svars = uIpms.Chan_base(ncs).vars{1};
            if num_svars > 0
                fignum = fignum + 1;
                figure(fignum)
                for nsvs = 1:num_svars
                    subplot(num_svars,numprot,numprot*(nsvs-1) + nprot)
                    plot(fitData.prot{nprot}.times_to_sim, outstruct.prot{nprot}.chan_memb_pred(ncs).sv{nsvs}(:,plotstep,plotseg))
                    title([uIpms.Chan_base(ncs).name ' membrane, state var ', uIpms.Chan_base(ncs).statenames{nsvs}, ', step ' num2str(plotstep) ', segment ' num2str(plotseg)])
                end
            end
        end
            
        
        if fitData.uIpms.fit.num_shells > 0

            fignum = 400;
            figure(fignum)
            subplot(1,numprot,nprot)
            plot(fitData.prot{nprot}.times_to_sim, outstruct.prot{nprot}.ImC(:,plotstep,plotseg), fitData.prot{nprot}.times_to_sim, squeeze(outstruct.prot{nprot}.ItC(:,plotstep,plotseg,:)))
            title(['Capacitive current, surface and each TT shell, step ' num2str(plotstep) ', segment ' num2str(plotseg)])
            
            for ncs = 1:num_chans
                fignum = fignum + 1;
                figure(fignum)
                subplot(1,numprot,nprot)
                plot(fitData.prot{nprot}.times_to_sim, outstruct.prot{nprot}.ImChan{ncs}(:,plotstep,plotseg), fitData.prot{nprot}.times_to_sim, squeeze(outstruct.prot{nprot}.ItChan{ncs}(:,plotstep,plotseg,:)))
                title([Chan_base(ncs).name ' current, surface and each TT shell, step ' num2str(plotstep) ', segment ' num2str(plotseg)])
            end
                

            fignum = fignum + 1;
            figure(fignum)
            subplot(1,numprot,nprot)
            net_TT_I = squeeze(outstruct.prot{nprot}.ItTot(:,plotstep,plotseg,:));
            net_memb_nonT = outstruct.prot{nprot}.ImTot - outstruct.prot{nprot}.ImT;
            plot(fitData.prot{nprot}.times_to_sim, net_memb_nonT(:,plotstep,plotseg), fitData.prot{nprot}.times_to_sim, net_TT_I)
            title(['Total current (cap plus channels), surface (nonT) and each TT shell, step ' num2str(plotstep) ', segment ' num2str(plotseg)])

            fignum = fignum + 1;
            figure(fignum)
            subplot(1,numprot,nprot)
            tot_TT_I = sum(net_TT_I,2);
            plot(fitData.prot{nprot}.times_to_sim, outstruct.prot{nprot}.ImT(:,plotstep,plotseg), '.', fitData.prot{nprot}.times_to_sim, tot_TT_I)
            title(['Total TT current, through Ra (dots) and through TT memb, step ' num2str(plotstep) ', segment ' num2str(plotseg)])

            
        end

        fignum = 499;
        for nms = 1:num_molecs
            if fitData.uIpms.fit.num_shells > 0
                fignum = fignum + 1;
                figure(fignum)
                subplot(1,numprot,nprot)
                plot(fitData.prot{nprot}.times_to_sim, squeeze(outstruct.prot{nprot}.molecs_t_rsh{nms}(:,plotstep,plotseg,:)))
                title([Molecs_base(nms).name ' tt concentration, each TT shell, step ' num2str(plotstep) ', segment ' num2str(plotseg)])
            end
            
            fignum = fignum + 1;
            figure(fignum)
            subplot(1,numprot,nprot)
            plot(fitData.prot{nprot}.times_to_sim, outstruct.prot{nprot}.molecs_i{nms}(:,plotstep,plotseg))
            title([Molecs_base(nms).name ' intrac concentration, step ' num2str(plotstep) ', segment ' num2str(plotseg)])
            
            fignum = fignum + 1;
            figure(fignum)
            subplot(1,numprot,nprot)
            plot(fitData.prot{nprot}.times_to_sim, outstruct.prot{nprot}.molecs_e{nms}(:,plotstep))
            title([Molecs_base(nms).name ' extrac concentration, step ' num2str(plotstep)])
            
        end
            
        fignum = fignum + 1;
        figure(fignum)
        if nprot == 1
            clf
        end
        subplot(1,numprot,nprot)
        hold on
        legendstr = [];
        concplot = 0;
        for nms = 1:num_molecs
            if fitData.uIpms.fit.num_shells > 0
                plot(fitData.prot{nprot}.times_to_sim, outstruct.prot{nprot}.molecs_t_rsh{nms}(:,plotstep,plotseg,plotshell))
                concplot = concplot + 1;
                legendstr{concplot} = [Molecs_base(nms).name 't'];
            end
            plot(fitData.prot{nprot}.times_to_sim, outstruct.prot{nprot}.molecs_i{nms}(:,plotstep,plotseg))
            concplot = concplot + 1;
            legendstr{concplot} = [Molecs_base(nms).name 'i'];
            plot(fitData.prot{nprot}.times_to_sim, outstruct.prot{nprot}.molecs_e{nms}(:,plotstep))
            concplot = concplot + 1;
            legendstr{concplot} = [Molecs_base(nms).name 'e'];
        end
        title(['all ion concentrations concentration, step ' num2str(plotstep)])
        legend(legendstr);    
            
        for nms = 1:num_molecs
            disp([Molecs_base(nms).name 'i = ' num2str(outstruct.prot{nprot}.molecs_i{nms}(end,plotstep,plotseg))]);
            disp([Molecs_base(nms).name 'e = ' num2str(outstruct.prot{nprot}.molecs_e{nms}(end,plotstep))]);
            if fitData.uIpms.fit.num_shells > 0
                disp([Molecs_base(nms).name 't = ' num2str(outstruct.prot{nprot}.molecs_t_rsh{nms}(end,plotstep,plotseg,plotshell))]);
            end
        end
        
        disp(['Ke*Cle = ' num2str(outstruct.prot{nprot}.molecs_e{2}(end,plotstep) * outstruct.prot{nprot}.molecs_e{3}(end,plotstep))]);
        disp(['Ki*Cli = ' num2str(outstruct.prot{nprot}.molecs_i{2}(end,plotstep) * outstruct.prot{nprot}.molecs_i{3}(end,plotstep))]);
        disp(['Nae*Ke*Cle^2 = ' num2str(outstruct.prot{nprot}.molecs_e{1}(end,plotstep) * outstruct.prot{nprot}.molecs_e{2}(end,plotstep) * outstruct.prot{nprot}.molecs_e{3}(end,plotstep).^2)]);
        disp(['Nai*Ki*Cli^2 = ' num2str(outstruct.prot{nprot}.molecs_i{1}(end,plotstep) * outstruct.prot{nprot}.molecs_i{2}(end,plotstep) * outstruct.prot{nprot}.molecs_i{3}(end,plotstep).^2)]);
        
        fignum = fignum + 1;
        figure(fignum)
        if nprot == 1
            clf
        end
        subplot(1,numprot,nprot)
        hold on
        legendstr = [];
        Nernst_ion_end = [];
        for nms = 1:num_molecs
            Nernst_ion = fitData.RTF.*log(outstruct.prot{nprot}.molecs_e{nms}(:,plotstep)./outstruct.prot{nprot}.molecs_i{nms}(:,plotstep,plotseg))./ Molecs_base(nms).valence;
            Nernst_ion_end(nms) = Nernst_ion(end);
            plot(fitData.prot{nprot}.times_to_sim, Nernst_ion)
            legendstr{nms} = Molecs_base(nms).name;
        end
        title('Nernst')
        legend(legendstr)
        hold off
        
        for nms = 1:num_molecs
            disp([Molecs_base(nms).name ' Nernst = ' num2str(Nernst_ion_end(nms))]);
        end
       
        for nms = 1:num_molecs
            fignum = fignum + 1;
            figure(fignum)
            if nprot == 1
                clf
            end
            subplot(1,numprot,nprot)
        
            total_intra_ion = sum(outstruct.prot{nprot}.molecs_i{nms}(:,plotstep,:) .* outstruct.prot{nprot}.volume_i_pred(:,plotstep,:),3);
            total_extra_ion = outstruct.prot{nprot}.molecs_e{nms}(:,plotstep) .* outstruct.prot{nprot}.volume_e_pred(:,plotstep);
            if uIpms.fit.num_shells > 0
                total_tt_ion = 0;
                for tkns = 1:uIpms.fit.num_shells
                    total_tt_ion = total_tt_ion + sum(outstruct.prot{nprot}.molecs_t_rsh{nms}(:,plotstep,:,tkns) .* outstruct.prot{nprot}.volume_t_pred(:,plotstep,:,tkns),3);
                end
                total_ion = total_intra_ion + total_extra_ion + total_tt_ion;
                plot(fitData.prot{nprot}.times_to_sim, total_intra_ion, fitData.prot{nprot}.times_to_sim, total_extra_ion, fitData.prot{nprot}.times_to_sim, total_tt_ion, fitData.prot{nprot}.times_to_sim, total_ion)
                legend('intra','extra','tt','total')
            else
                total_ion = total_intra_ion + total_extra_ion;
                plot(fitData.prot{nprot}.times_to_sim, total_intra_ion, fitData.prot{nprot}.times_to_sim, total_extra_ion, fitData.prot{nprot}.times_to_sim, total_ion)
                legend('intra','extra','total')
            end
            title(['total ' Molecs_base(nms).name ' amounts'])
        end
        
        total_cap = fitData.tot_surf_area .* all_params(1);
        if fitData.uIpms.fit.num_shells > 0
            disp(['current due to voltage across Ra (excludes diffusion current) = ' num2str(outstruct.prot{nprot}.Itot_m_T(end,plotstep)./ total_cap)]);
        end
        disp(['membrane cap current = ' num2str(outstruct.prot{nprot}.Itot_m_C(end,plotstep)./ total_cap)]);
        for ncs = 1:num_chans
            disp(['membrane ' Chan_base(ncs).name ' current = ' num2str(outstruct.prot{nprot}.Itot_m_chan{ncs}(end,plotstep)./ total_cap)]);
        end
        
        if fitData.uIpms.fit.num_shells > 0
            for ncs = 1:num_chans
                disp(['ttubule ' Chan_base(ncs).name ' current = ' num2str(outstruct.prot{nprot}.Itot_t_chan{ncs}(end,plotstep)./ total_cap)]);
            end
        end

        fignum = 599;
        if fitData.uIpms.fit.num_shells > 0
            fignum = fignum + 1;
            figure(fignum)
            if nprot == 1
                clf
            end
            subplot(1,numprot,nprot)
            plot(fitData.prot{nprot}.times_to_sim,outstruct.prot{nprot}.volume_t_tot(:,plotstep))
            title('total t-tubule volume')
        end

        fignum = fignum + 1;
        figure(fignum)
        if nprot == 1
            clf
        end
        subplot(1,numprot,nprot)
        plot(fitData.prot{nprot}.times_to_sim,outstruct.prot{nprot}.volume_i_tot(:,plotstep))
        title('total intracellular volume')

        fignum = fignum + 1;
        figure(fignum)
        if nprot == 1
            clf
        end
        subplot(1,numprot,nprot)
        plot(fitData.prot{nprot}.times_to_sim,outstruct.prot{nprot}.volume_e_pred(:,plotstep))
        title('total extracellular volume')

        
        if nprot == 1
            fignum = 699;
            
            V = [-120:60];
            
            vouts_chan = [];
            for ncs = 1:num_chans
                tc = uIpms.Chan_base(ncs);
                num_v_outs = tc.vars{4};
                if num_v_outs > 0
                    fignum = fignum + 1;
                    figure(fignum)
                    mol_array = [];
                    for nms = 1:length(tc.molecs_needed)
                        molnum = tc.molecs_needed(nms);
                        mol_array{nms,1} = outstruct.prot{nprot}.molecs_i{molnum}(end,plotstep,plotseg);
                        mol_array{nms,2} = outstruct.prot{nprot}.molecs_e{molnum}(end,plotstep);
                    end
                    z = Molecs_base(tc.ion_current_mult(1,1)).valence;
                    
                    if uIpms.each_step_diff_exp == 1
                        numprmshere = length(uIpms.fit.param_list_idxnum{ncs}{nprot});
                        for prmnum = 1:numprmshere
                            uIpms.chan_param_idxs{ncs}(prmnum) = uIpms.fit.param_list_idxnum{ncs}{nprot}{prmnum}{plotstep};
                        end
                    end
                    
                    [vouts_chan{ncs}.vo, vout_names] = tc.ChanFunc([],V,all_params(uIpms.chan_param_idxs{ncs}),4,[],[],fitData.RTF, mol_array, z, tc.extra_params); 
                    for nvo = 1:num_v_outs
                        subplot(ceil(num_v_outs/2),2,nvo)
                        plot(V,vouts_chan{ncs}.vo{nvo})
                        title([tc.name ' ' vout_names{nvo}])
                    end
                end
            end
                
            % I vs V at steady state for each channel
            fignum = 799;
            fignum = fignum + 1;
            figure(fignum)
            if nprot == 1
                clf
            end
            subplot(1,numprot,nprot)
            
            hold on
            legendstr = [];

            for ncs = 1:num_chans
                tc = uIpms.Chan_base(ncs);
                ss_states_idxs = tc.vars{5};
                chan_states = [];
                for ssnum = 1:length(ss_states_idxs)
                    chan_states = [chan_states ,vouts_chan{ncs}.vo(ss_states_idxs(ssnum))];
                end
                thischan_params = all_params(uIpms.chan_param_idxs{ncs});
                thischan_perm = thischan_params(1);
                thischan_partial_cond_factor = thischan_perm .* fitData.sarc_seg_area .* 1000;
                mol_array = [];
                for nms = 1:length(tc.molecs_needed)
                    molnum = tc.molecs_needed(nms);
                    mol_array{nms,1} = outstruct.prot{nprot}.molecs_i{molnum}(end,plotstep,plotseg);
                    mol_array{nms,2} = outstruct.prot{nprot}.molecs_e{molnum}(end,plotstep);
                end
                z = Molecs_base(tc.ion_current_mult(1,1)).valence;

                thischan_max_cond_factor = find_max_cond_factor(tc.conductance_constant, thischan_partial_cond_factor, mol_array{1,1}, mol_array{1,2}, z, tc.use_GHK, fitData.F, fitData.RTF);
                [~, ImoChan{ncs}] = tc.ChanFunc(chan_states, V, all_params(uIpms.chan_param_idxs{ncs}), 2, thischan_max_cond_factor, tc.use_GHK, fitData.RTF, mol_array, z, tc.extra_params);
                plot(V,ImoChan{ncs})
                legendstr{ncs} = tc.name;
            end
            
            title(['Steady state currents vs. V, sarc membrane only'])
            legend(legendstr)
            
            fignum = fignum + 1;
            figure(fignum)
            Isum_sarc = 0;
            for ncs = 1:num_chans
                Isum_sarc = Isum_sarc + ImoChan{ncs};
            end
            plot(V, Isum_sarc)
            title('resting currents, sum of all channels')
            
            % plot selected currents at final time point
            fignum = fignum + 1;
            figure(fignum)
            if nprot == 1
                clf
            end
            subplot(1,numprot,nprot)
            hold on
            legendstr = [];
            currsel = [1 2 3 4 5 9];
            numc = 0;
            for cs = currsel
                numc = numc + 1;
                tc = Chan_base(cs);
                plot(V,ImoChan{cs})
                legendstr{numc} = tc.name;
            end

            % add zero current
            tempzero = zeros(1,length(V));
            plot(V,tempzero,'k')

            xline(outstruct.prot{nprot}.Vm_pred(end,i,fitData.uIpms.gen.voltage_probe_segment))

            title(['Final steady state currents vs. V, sarc membrane only'])
            legend(legendstr)


            % plot selected currents at specified time point

            chosen_time_index = 29999;
            %chosen_time_index = 30011;
            %chosen_time_index = 59000;
            %chosen_time_index = 33000;
            %chosen_time_index = 30500;
            %chosen_time_index = 30500;
            for cs = 1:num_chans
                tc = uIpms.Chan_base(cs);

                num_v_outs = tc.vars{4};
                if num_v_outs > 0
                    mol_array = [];
                    for nms = 1:length(tc.molecs_needed)
                        molnum = tc.molecs_needed(nms);
                        mol_array{nms,1} = outstruct.prot{nprot}.molecs_i{molnum}(chosen_time_index,plotstep,plotseg);
                        mol_array{nms,2} = outstruct.prot{nprot}.molecs_e{molnum}(chosen_time_index,plotstep);
                    end
                    z = Molecs_base(tc.ion_current_mult(1,1)).valence;
                    
                    if uIpms.each_step_diff_exp == 1
                        numprmshere = length(uIpms.fit.param_list_idxnum{cs}{nprot});
                        for prmnum = 1:numprmshere
                            uIpms.chan_param_idxs{cs}(prmnum) = uIpms.fit.param_list_idxnum{cs}{nprot}{prmnum}{plotstep};
                        end
                    end
                    
                    [vouts_chan{cs}.vo, vout_names] = tc.ChanFunc([],V,all_params(uIpms.chan_param_idxs{cs}),4,[],[],fitData.RTF, mol_array, z, tc.extra_params); 

                else

                    mol_array = [];
                    for nms = 1:length(tc.molecs_needed)
                        molnum = tc.molecs_needed(nms);
                        mol_array{nms,1} = outstruct.prot{nprot}.molecs_i{molnum}(chosen_time_index,plotstep,plotseg);
                        mol_array{nms,2} = outstruct.prot{nprot}.molecs_e{molnum}(chosen_time_index,plotstep);
                    end
                    z = Molecs_base(tc.ion_current_mult(1,1)).valence;

                end

                ss_states_idxs = tc.vars{5};
                chan_states = [];
                for ssnum = 1:length(ss_states_idxs)
                    chan_states = [chan_states ,vouts_chan{cs}.vo(ss_states_idxs(ssnum))];
                end
                thischan_params = all_params(uIpms.chan_param_idxs{cs});
                thischan_perm = thischan_params(1);
                thischan_partial_cond_factor = thischan_perm .* fitData.sarc_seg_area .* 1000;
                    
                thischan_max_cond_factor = find_max_cond_factor(tc.conductance_constant, thischan_partial_cond_factor, mol_array{1,1}, mol_array{1,2}, z, tc.use_GHK, fitData.F, fitData.RTF);
                [~, ImoChan{cs}] = tc.ChanFunc(chan_states, V, all_params(uIpms.chan_param_idxs{cs}), 2, thischan_max_cond_factor, tc.use_GHK, fitData.RTF, mol_array, z, tc.extra_params);

            end

            fignum = fignum + 1;
            figure(fignum)
            if nprot == 1
                clf
            end
            subplot(1,numprot,nprot)
            hold on

            legendstr = [];
            %currsel = [1 2 3 4 5 9 10];
            currsel = [5 9];
            numc = 0;

            for cs = currsel
                numc = numc + 1;
                tc = uIpms.Chan_base(cs);

                plot(V,ImoChan{cs})
                legendstr{numc} = tc.name;
            end

            % add zero current
            tempzero = zeros(1,length(V));
            plot(V,tempzero,'k')

            xline(outstruct.prot{nprot}.Vm_pred(chosen_time_index,i,fitData.uIpms.gen.voltage_probe_segment))

            title(['Steady state currents vs. V at time index, sarc membrane only'])
            legend(legendstr)


            tot_current_allmols = 0;
            for nms = 1:num_molecs
                
                Molrate_current_vs_V{nms} = 0;
    
                for nchs = 1:size(uIpms.molec_chan{nms},1)
                    thischan = uIpms.molec_chan{nms}(nchs,1);
                    thismult = uIpms.molec_chan{nms}(nchs,2);
                    Molrate_current_vs_V{nms} = Molrate_current_vs_V{nms} + thismult.*ImoChan{thischan};
                    if fitData.uIpms.fit.num_shells > 0
                        Molrate_current_vs_V{nms} = Molrate_current_vs_V{nms} + thismult.*outstruct.prot{nprot}.Itot_t_chan{thischan}(:,plotstep);
                    end
                end

                tot_current_allmols = tot_current_allmols + Molrate_current_vs_V{nms};
            end
                
            fignum = fignum + 1;
            figure(fignum)
            if nprot == 1
                clf
            end
            subplot(1,numprot,nprot)
            hold on
            legendstr = [];

            plot(V, tot_current_allmols)
            legendstr{1} = 'total';

            for nms = 1:num_molecs
                plot(V, Molrate_current_vs_V{nms})
                legendstr{nms+1} = Molecs_base(nms).name;
            end

            xline(outstruct.prot{nprot}.Vm_pred(chosen_time_index,i,fitData.uIpms.gen.voltage_probe_segment))
            title(['Net total or ion current across surface membrane vs V, at specified time'])
            legend(legendstr)
        

            if 1

                chosen_time_indices = [29999 30050 59999];

                for cti = chosen_time_indices

                    disp([' '])
                    disp(['time index ' num2str(cti)])

                    disp(['Vm = ' num2str(outstruct.prot{nprot}.Vm_pred(cti,1,fitData.uIpms.gen.voltage_probe_segment))])
                
                    for nms = 1:num_molecs
                        disp([Molecs_base(nms).name 'i = ' num2str(outstruct.prot{nprot}.molecs_i{nms}(cti,plotstep,plotseg))]);
                        disp([Molecs_base(nms).name 'e = ' num2str(outstruct.prot{nprot}.molecs_e{nms}(cti,plotstep))]);
                        if fitData.uIpms.fit.num_shells > 0
                            disp([Molecs_base(nms).name 't = ' num2str(outstruct.prot{nprot}.molecs_t_rsh{nms}(cti,plotstep,plotseg,plotshell))]);
                        end
                    end
                    
                    for nms = 1:num_molecs
                        Nernst_ion_vs_t = fitData.RTF.*log(outstruct.prot{nprot}.molecs_e{nms}(:,plotstep)./outstruct.prot{nprot}.molecs_i{nms}(:,plotstep,plotseg))./ Molecs_base(nms).valence;
                        disp([Molecs_base(nms).name ' Nernst = ' num2str(Nernst_ion_vs_t(cti))]);
                    end

                    for nms = 1:num_molecs
                        disp(['Net ' Molecs_base(nms).name 'current (all channels) ' num2str(Molrate_current{nms}(cti)/total_cap)])
                    end
                    
                end

            end


        end




        
        fignum = 899;
        if isfield(fitData.prot{nprot},'do_obj_func_calcs') && fitData.prot{nprot}.do_obj_func_calcs == 1 
            
            disp(['Value of Objective function = ' num2str(outstruct.obj_val.net)]);

            fignum = fignum + 1;
            figure(fignum)
            plot(outstruct.obj_residuals.prot{nprot}.resid)
            title('residuals')
            
            if isequal(protocol{nprot}.obj_func, @AP_phase_plane)
                fignum = fignum + 1;
                figure(fignum)
                if nprot == 1
                    clf
                end
                subplot(1,numprot,nprot)

                hold off
                plot(outstruct.extra_obj_out.prot{nprot}.V_binned_pos, outstruct.extra_obj_out.prot{nprot}.dVdt_binned_pos_expmt,'b.')
                hold on
                plot(outstruct.extra_obj_out.prot{nprot}.V_binned_pos, outstruct.extra_obj_out.prot{nprot}.dVdt_binned_pos_pred,'r.')
                plot(outstruct.extra_obj_out.prot{nprot}.V_binned_neg, outstruct.extra_obj_out.prot{nprot}.dVdt_binned_neg_expmt,'b.')
                plot(outstruct.extra_obj_out.prot{nprot}.V_binned_neg, outstruct.extra_obj_out.prot{nprot}.dVdt_binned_neg_pred,'r.')
                
            end
            
        end
        
        

    end     % if do_medium_graphs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end     % for nprot = 1:numprot

        

