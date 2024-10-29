function chan_perm_factors = set_chan_perm_factors(protocol, fitData)

chan_perm_factors = [];

for ncs = 1:fitData.uIpms.num_chans
    chan_perm_factors(ncs) = 1;
    if isfield(protocol, 'modified_perm_chans')
        chanlist = protocol.modified_perm_chans;
        idx = find(chanlist == fitData.uIpms.Chan_base(ncs).name,1);
        if ~isempty(idx)
            chan_perm_factors(ncs) = protocol.modified_perms(idx);
        end
    end
end
        
    
    