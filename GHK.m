 function GHK_out = GHK(Vm, co, ci, z, Vprime, RTF)
 % Vprime is to account for surface charge layers that alter conductance compared to that expected due to bulk Nernst
 
 % GHK_out is in units of mV*mM, a measure of the driving force for current flow.
 
if Vm == Vprime
    GHK_out = RTF * (ci - co .* exp(-z.*Vprime./RTF) .* exp(-z.*(Vm - Vprime)./RTF));
else
    GHK_out = (Vm - Vprime) .* (ci - co .* exp(-z.*Vm./RTF)) ./ (1 - exp(-z.*(Vm - Vprime)./RTF));
end


