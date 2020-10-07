%%
function lfit_val = fit_output(temp,w_0,p,nu_eff,w_r)
% calculate the fit output values given the parameters and temperatures
    k_b = 0.695;
    lfit_val = w_0.*(1./(1-exp(-nu_eff./(k_b.*temp)))).^p + w_r;
    
end