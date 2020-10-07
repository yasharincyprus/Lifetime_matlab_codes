%%
function lfit_val = fit_output(temp,w_0,w_1,w_2,p,q,r,nu_eff,nu_eff1,nu_eff2,w_r)
% calculate the fit output values given the parameters and temperatures
    k_b = 0.695;
    lfit_val = w_0.*(1+((1./(-1+exp(nu_eff./(k_b.*temp))))).^p)+w_1.*(1+((1./(-1+exp(nu_eff1./(k_b.*temp))))).^q)+w_2.*(1+((1./(-1+exp(nu_eff2./(k_b.*temp))))).^r) + w_r;
    
end