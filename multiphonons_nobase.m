clearvars; clc; close all;
data_mat = csvread('transrate_temp_singexp_nofilt.csv');
temperature = data_mat(:,1);
w_t = data_mat(:,2);

E_g = 5529;
w_r = 0;

p = 4;
q = 3;
r = 1;

% nu_eff = floor(E_g / p);

par(1,1) = 100;
par(2,1) = 100;
par(3,1) = 3000;
% par(4,1) = 0;


% for p = 1:20
% nu_eff = floor(E_g / p);
nu_eff = 884;
nu_eff1 = 595;
nu_eff2 = 128;

tic;
options=optimset('Display','final','TolFun',1e-9,'TolX',1e-10,...
                'MaxFunEvals',5e4,'MaxIter',1e6);
            
Niter=1;maxiter=50;ssq=100*ones(length(maxiter),1);ep1=1e-5;

while Niter<50 && ssq(Niter)>0.001
    Niter=Niter+1;
   [parmin,fval,exitflag]=fminsearch(@ssqmin_nobase,par,options,w_t,temperature,p,q,r,nu_eff,nu_eff1,nu_eff2,w_r);
    npar=parmin;
    fit=fit_output_nobase(temperature,npar(1,1),npar(2,1),npar(3,1),p,q,r,nu_eff,nu_eff1,nu_eff2,w_r);
    ssq(Niter)=sum((fit-w_t).^2);
    if (ssq(Niter)-ssq(Niter-1))/ssq(Niter-1)>-ep1
        break
    end
    t1=toc;
    disp(['iteration ',num2str(Niter-1),',time ',num2str(t1),' secs,ssq ',num2str(ssq(Niter))])
end

figure
scatter(temperature,w_t,'o','k')
hold on
plot(temperature,fit,'-k')
% end
