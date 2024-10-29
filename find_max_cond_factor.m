function max_cond_factor = find_max_cond_factor(channel_constant, partial_factor, ci, co, z, use_GHK, F, RTF)

if isempty(channel_constant)
    if use_GHK
        ion_factor = F .* z.^2 ./ RTF;
    else
        Nernst = RTF.*log(co./ci)./z;
        ion_factor = F .* Nernst .* co .* z.^3 ./ (RTF.^2 .* (exp(z.*Nernst./RTF) - 1));
    end
else
    ion_factor = channel_constant; 
end

max_cond_factor = ion_factor .* partial_factor;

