if fitData.uIpms.fit.use_IC_comps && (sP.num_segs >= 2 || sP.num_shells >= 2)

    gir = sP.gIR;
    gil = sP.gIL;
    gtr = sP.gT;
    gtl = sP.gTL;
    gs = sP.g_S;
    ct = sP.C_t;
    ra = sP.ra;

    if ~sP.try_small_matrix
        sP.A = zeros(2*sP.num_segs*sP.num_shells, 2*sP.num_segs*sP.num_shells);
        sP.B = zeros(sP.num_segs*sP.num_shells, sP.num_segs*sP.num_shells);

        for n = 1:sP.num_shells
            if n == 1
                nm1 = 0;
            else
                nm1 = 1;
            end
            if n == sP.num_shells
                nmn = 0;
                nmf = 1;
            else
                nmn = 1;
                nmf = 0;
            end
            for j = 1:sP.num_segs

                rw = (n-1)*sP.num_segs + j;

                if j == 1
                    jm1 = 0;
                else
                    jm1 = 1;
                end
                if j == sP.num_segs
                    jmj = 0;
                else
                    jmj = 1;
                end
                if j == 1 && j == sP.num_segs
                    jmf = 0;
                elseif j == 1 || j == sP.num_segs
                    jmf = 1;
                else
                    jmf = 2;
                end


                if nm1
                    sP.A(rw,(n-2)*ns+j) = gir(n-1);
                end
                if jm1
                    sP.A(rw,(n-1)*ns+j-1) = gil(n);
                end

                sP.A(rw,(n-1)*ns+j) = - jmf*gil(n) - nmf*sP.g_S;
                if nm1
                    sP.A(rw,(n-1)*ns+j) = sP.A(rw,(n-1)*ns+j) - gir(n-1);
                end
                if nmn
                    sP.A(rw,(n-1)*ns+j) = sP.A(rw,(n-1)*ns+j) - gir(n);
                end

                if jmj
                    sP.A(rw,(n-1)*ns+j+1) = gil(n);
                end
                if nmn
                    sP.A(rw,n*ns+j) = gir(n);
                end


                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if nm1
                    sP.B(rw,(n-2)*ns+j) = gtr(n-1);
                end
                if jm1
                    sP.B(rw,(n-1)*ns+j-1) = gtl(n);
                end

                sP.B(rw,(n-1)*ns+j) = -jmf*gtl(n) - nmf/ra;
                if nm1
                    sP.B(rw,(n-1)*ns+j) = sP.B(rw,(n-1)*ns+j) - gtr(n-1);
                end
                if nmn
                    sP.B(rw,(n-1)*ns+j) = sP.B(rw,(n-1)*ns+j) - gtr(n);
                end

                if jmj
                    sP.B(rw,(n-1)*ns+j+1) = gtl(n);
                end
                if nmn
                    sP.B(rw,n*ns+j) = gtr(n);
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                rw = rw + ns*nsh;

                if nm1
                    sP.A(rw,(n-2)*ns+j) = gtr(n-1);
                end
                if jm1
                    sP.A(rw,(n-1)*ns+j-1) = gtl(n);
                end

                sP.A(rw,(n-1)*ns+j) = - jmf*gtl(n) - nmf/ra;
                if nm1
                    sP.A(rw,(n-1)*ns+j) = sP.A(rw,(n-1)*ns+j) - gtr(n-1);
                end
                if nmn
                    sP.A(rw,(n-1)*ns+j) = sP.A(rw,(n-1)*ns+j) - gtr(n);
                end

                if jmj
                    sP.A(rw,(n-1)*ns+j+1) = gtl(n);
                end
                if nmn
                    sP.A(rw,n*ns+j) = gtr(n);
                end


            end

        end

        sP.A(1:ns*nsh, ns*nsh+1:2*ns*nsh) = diag(-ct);
        sP.A(ns*nsh+1:2*ns*nsh, ns*nsh+1:2*ns*nsh) = diag(ct);

    end


    if sP.try_small_matrix

        sP.A = zeros(sP.num_segs*sP.num_shells, sP.num_segs*sP.num_shells);
        sP.B = zeros(sP.num_segs*sP.num_shells, sP.num_segs*sP.num_shells);
        sP.TT = zeros(sP.num_segs*sP.num_shells, sP.num_segs*sP.num_shells);

        for n = 1:sP.num_shells
            if n == 1
                nm1 = 0;
            else
                nm1 = 1;
            end
            if n == sP.num_shells
                nmn = 0;
                nmf = 1;
            else
                nmn = 1;
                nmf = 0;
            end
            for j = 1:sP.num_segs

                rw = (n-1)*sP.num_segs + j;

                if j == 1
                    jm1 = 0;
                else
                    jm1 = 1;
                end
                if j == sP.num_segs
                    jmj = 0;
                else
                    jmj = 1;
                end
                if j == 1 && j == sP.num_segs
                    jmf = 0;
                elseif j == 1 || j == sP.num_segs
                    jmf = 1;
                else
                    jmf = 2;
                end


                if nm1
                    sP.A(rw,(n-2)*ns+j) = sP.gIR(n-1) + sP.gT(n-1);
                end
                if jm1
                    sP.A(rw,(n-1)*ns+j-1) = sP.gIL(n) + sP.gTL(n);
                end

                sP.A(rw,(n-1)*ns+j) = - jmf*sP.gIL(n) - jmf*sP.gTL(n) - nmf*sP.g_S - nmf/sP.ra;
                if nm1
                    sP.A(rw,(n-1)*ns+j) = sP.A(rw,(n-1)*ns+j) - sP.gIR(n-1) - sP.gT(n-1);
                end
                if nmn
                    sP.A(rw,(n-1)*ns+j) = sP.A(rw,(n-1)*ns+j) - sP.gIR(n) - sP.gT(n);
                end

                if jmj
                    sP.A(rw,(n-1)*ns+j+1) = sP.gIL(n) + sP.gTL(n);
                end
                if nmn
                    sP.A(rw,n*ns+j) = sP.gIR(n) + sP.gT(n);
                end


                if nm1
                    sP.B(rw,(n-2)*ns+j) = sP.gT(n-1);
                end
                if jm1
                    sP.B(rw,(n-1)*ns+j-1) = sP.gTL(n);
                end

                sP.B(rw,(n-1)*ns+j) = -jmf*sP.gTL(n) - nmf/sP.ra;
                if nm1
                    sP.B(rw,(n-1)*ns+j) = sP.B(rw,(n-1)*ns+j) - sP.gT(n-1);
                end
                if nmn
                    sP.B(rw,(n-1)*ns+j) = sP.B(rw,(n-1)*ns+j) - sP.gT(n);
                end

                if jmj
                    sP.B(rw,(n-1)*ns+j+1) = sP.gTL(n);
                end
                if nmn
                    sP.B(rw,n*ns+j) = sP.gT(n);
                end



                if nm1
                    sP.TT(rw,(n-2)*ns+j) = -sP.gT(n-1);
                end
                if jm1
                    sP.TT(rw,(n-1)*ns+j-1) = -sP.gTL(n);
                end

                sP.TT(rw,(n-1)*ns+j) = jmf*sP.gTL(n) + nmf/sP.ra;
                if nm1
                    sP.TT(rw,(n-1)*ns+j) = sP.TT(rw,(n-1)*ns+j) + sP.gT(n-1);
                end
                if nmn
                    sP.TT(rw,(n-1)*ns+j) = sP.TT(rw,(n-1)*ns+j) + sP.gT(n);
                end

                if jmj
                    sP.TT(rw,(n-1)*ns+j+1) = -sP.gTL(n);
                end
                if nmn
                    sP.TT(rw,n*ns+j) = -sP.gT(n);
                end


            end
        end
    end


    if sP.input_method == 0 || sP.input_method == 4 || sP.input_method == 5 || sP.input_method == 6
        if ~sP.try_small_matrix
            sP.d_eqn_corr = zeros(2*sP.num_segs*sP.num_shells,1);
        else
            sP.d_eqn_corr = zeros(sP.num_segs*sP.num_shells,1);
        end
        vpssh_idx = (sP.voltage_probe_shell - 1).* ns + sP.voltage_probe_seg;
        nonzero_rows = find(sP.A(:,vpssh_idx));
        for nzr = nonzero_rows
            sP.d_eqn_corr(nzr) = sP.A(nzr,vpssh_idx);
        end

        sP.A(vpssh_idx,:) = [];
        sP.A(:,vpssh_idx) = [];
    end

end

if ~fitData.uIpms.fit.use_IC_comps && sP.num_segs >= 2

    sP.B = sP.gIL .* (diag(2.*ones(1,sP.num_segs)) + diag(-1.*ones(1,sP.num_segs-1),+1) + diag(-1.*ones(1,sP.num_segs-1),-1) );
    sP.B(1,1) = sP.B(1,1)/2;
    sP.B(end,end) = sP.B(end,end)/2;

    sP.A = -sP.B .* sP.C ./ sP.g_S;
    if nsh > 0
        sP.A = sP.A - diag( (sP.C + sP.C/(sP.g_S*sP.ra)) .* ones(1,sP.num_segs));
    else
        sP.A = sP.A - diag(sP.C .* ones(1,sP.num_segs));
    end

    if sP.input_method == 0 || sP.input_method == 4 || sP.input_method == 5 || sP.input_method == 6
        sP.d_eqn_corr_above = sP.A(sP.voltage_probe_seg-1,sP.voltage_probe_seg);    % for when 1 equation is removed when using V input
        sP.d_eqn_corr_below = sP.A(sP.voltage_probe_seg+1,sP.voltage_probe_seg);

        sP.A(sP.voltage_probe_seg,:) = [];
        sP.A(:,sP.voltage_probe_seg) = [];
    end


end
