function [Vc, dVc_dt]  = sigmoid_step_func(t,rise,offset,magn,init)

toobig = rise.*(t-offset);

if any(toobig > 100,'all')
    t_vec = t(:);
    toobig_size = size(toobig);
    toobig_vec = toobig(:);
    numtb = length(toobig_vec);
    keyval_vec = ones(numtb,1);
    for tb = 1:numtb
        if toobig(tb) > 100
            keyval_vec(tb) = 1;
        else
            keyval_vec(tb) = exp(rise.*(t_vec(tb)-offset)) ./ (1 + exp(rise.*(t_vec(tb)-offset)));
        end
    end
    keyval = reshape(keyval_vec,toobig_size);
else
    keyval = exp(rise.*(t-offset)) ./ (1 + exp(rise.*(t-offset)));
end

Vc = (magn-init).*keyval+init;

dVc_dt = (magn-init).*rise.*keyval.*(1-keyval);


