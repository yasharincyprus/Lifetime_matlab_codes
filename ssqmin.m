%%
function ssq = peakfitmin(par,rate,temp,p,nu_eff,w_r)
% controls how ssq is minimised. 

    lfit_val = fit_output(temp,par(1,1),p,nu_eff,w_r);
    ssq = sum((lfit_val - rate).^2);
    
    % constrains absolute zero rate to be positive and smaller than 10-K
    % rate
    
    pgf = 1e20*(par(1,1)<0);
    glf = 1e20*(par(1,1)>min(rate));
    ssq = ssq*(1 + (par(1,1)^2)*(pgf+glf));
    
%     pgz = 1e20*(par(2,1)<0);
%     glz = 1e20*(par(2,1)>min(rate));
%     ssq = ssq*(1 + (par(2,1)^2)*(pgz+glz));
    
end
