classdef MolecClass
    properties
        valence = 1;
        diffusion_constant = 1e-5;  % cm^2/s
        init_tt = 100;  % initial t-tubule concentration, mM
        init_intrac = 10;   % initial intracellular concentration, mM
        init_extrac = 100;  % initial extracellular concentration, mM
        name
    end
    
    methods
        function V = Nernst(obj,ce,ci,RTF)
            V = RTF.*log(ce./ci)./obj.valence;
        end
    end
end