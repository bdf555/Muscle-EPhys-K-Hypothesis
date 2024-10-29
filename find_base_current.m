function I = find_base_current(Vm, cond_factor, co, ci, use_GHK, Vprime, z, RTF)
% cond_factor is either conductance (uS) if use_GHK = 0, or conductance per mM (uS/mM) if use_GHK = 1
% co and ci are extracellular and intracellular concentration in mM
% z is valence of ion

% I in nA

if use_GHK
    
    ghk = GHK(Vm, co, ci, z, Vprime, RTF);
    I = cond_factor .* ghk;
    
else
    
    Nernst = RTF.*log(co./ci)./z;
    I = (Vm - Nernst).*cond_factor;
    
end