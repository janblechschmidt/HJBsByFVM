clear all
% tau = 0.5;
H = 40;
V = 40;
F = 40;
R = 0.05;
sigmaH = 0.3;
sigmaV=0.3;
rho_VH = 0.5;

Tmax = 1.0;
gamma1 = @(V,H,F,tau) (log(H/F) + (R - .5*sigmaH^2)*tau) / (sigmaH*sqrt(tau));
gamma2 = @(V,H,F,tau) (log(V/F) + (R - .5*sigmaV^2)*tau) / (sigmaV*sqrt(tau));
sigma = sqrt( sigmaV^2 + sigmaH^2 - 2*rho_VH*sigmaV*sigmaH );

alpha1 = @(V,H,F,tau) gamma1(V,H,F,tau)+sigmaH*sqrt(tau);
alpha2 = @(V,H,F,tau) ( log(V/H) - 0.5*sigma^2*tau ) / (sigma * sqrt(tau));
% alpha2 = @(V,H,F,tau) ( log(V/H) - 0.5*sigma^2*sqrt(tau) ) / (sigma * sqrt(tau));
rho_H = (rho_VH*sigmaV-sigmaH)/sigma;

beta1 = @(V,H,F,tau) gamma2(V,H,F,tau)+sigmaV*sqrt(tau);
beta2 = @(V,H,F,tau) ( log(H/V) - 0.5*sigma^2*tau ) / (sigma * sqrt(tau));
% beta2 = @(V,H,F,tau) ( log(H/V) - 0.5*sigma^2*sqrt(tau) ) / (sigma * sqrt(tau));
rho_V = (rho_VH*sigmaH-sigmaV)/sigma;

% M = @(x,F,tau) x(:,2).*mvncdf([alpha1, alpha2],0,rho_H) ...
%     + x(:,1).*mvncdf([beta1,beta2],0,rho_V) ...
%     - F*exp(-R*tau).*mvncdf([gamma1,gamma2],0,rho_VH);
% covmat = @(rho) [1, -sqrt(1-rho^2);  -sqrt(1-rho^2), 1];
covmat = @(rho) [1, rho;  rho, 1];
M = @(V,H,F,tau) H.*mvncdf([alpha1(V,H,F,tau), alpha2(V,H,F,tau)],[0,0],covmat(rho_H)) ...
    + V.*mvncdf([beta1(V,H,F,tau),beta2(V,H,F,tau)],[0,0],covmat(rho_V)) ...
    - F*exp(-R*tau).*mvncdf([gamma1(V,H,F,tau),gamma2(V,H,F,tau)],[0,0],covmat(rho_VH));

solfunc = @(t,x) exp(-R*(Tmax-t)) * F - M(x(:,1),x(:,2),0,Tmax-t) + M(x(:,1),x(:,2),F,Tmax-t);
solfunc(0.5,[40,40])
% 
% M2 = @(V,H,F,tau) H.*bivnormcdf(alpha1(V,H,F,tau), alpha2(V,H,F,tau),rho_H) ...
%     + V.*bivnormcdf(beta1(V,H,F,tau),beta2(V,H,F,tau),rho_V) ...
%     - F*exp(-R*tau).*bivnormcdf(gamma1(V,H,F,tau),gamma2(V,H,F,tau),rho_VH);
% 
% 
% solfunc2 = @(t,x) exp(-R*(Tmax-t)) * F - M2(x(:,1),x(:,2),0.0001,Tmax-t) + M2(x(:,1),x(:,2),F,Tmax-t);
% solfunc2(0.5,[40,40])

% TODO: Plot function
[X,Y]=meshgrid(linspace(0,100,101), linspace(0,100,101));

V=40.0;H=40.0001;F=0;tau=1e-8;
M(V,H,F,tau)
V=40.0;H=0.0000000001;F=40;tau=0.6;
M(V,H,F,tau)
solfunc(Tmax-tau,[V,H])*exp(R*tau)
