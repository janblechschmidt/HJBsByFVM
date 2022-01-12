function solfunc = explicitSolutionWorstOfTwoAssetsPut(F,R,sigmaH,sigmaV,rho_VH,Tmax)
myLog = @(x,y) max(log(x./y),-inf);
gamma1 = @(V,H,F,tau) (myLog(H,F) + (R - .5*sigmaH^2)*tau) / (sigmaH*sqrt(tau));
gamma2 = @(V,H,F,tau) (myLog(V,F) + (R - .5*sigmaV^2)*tau) / (sigmaV*sqrt(tau));
sigma = sqrt( sigmaV^2 + sigmaH^2 - 2*rho_VH*sigmaV*sigmaH );

alpha1 = @(V,H,F,tau) gamma1(V,H,F,tau)+sigmaH*sqrt(tau);
alpha2 = @(V,H,F,tau) ( myLog(V,H) - 0.5*sigma^2*tau ) / (sigma * sqrt(tau));
rho_H = (rho_VH*sigmaV-sigmaH)/sigma;

beta1 = @(V,H,F,tau) gamma2(V,H,F,tau)+sigmaV*sqrt(tau);
beta2 = @(V,H,F,tau) ( myLog(H,V) - 0.5*sigma^2*tau ) / (sigma * sqrt(tau));
rho_V = (rho_VH*sigmaH-sigmaV)/sigma;

covmat = @(rho) [1, rho;  rho, 1];
M = @(V,H,F,tau) H.*mvncdf([alpha1(V,H,F,tau), alpha2(V,H,F,tau)],0,covmat(rho_H)) ...
    + V.*mvncdf([beta1(V,H,F,tau),beta2(V,H,F,tau)],0,covmat(rho_V)) ...
    - F*exp(-R*tau).*mvncdf([gamma1(V,H,F,tau),gamma2(V,H,F,tau)],0,covmat(rho_VH));

solfunc = @(t,x) exp(-R*(Tmax-t)) * F - M(x(:,1),x(:,2),0,Tmax-t) + M(x(:,1),x(:,2),F,Tmax-t);
