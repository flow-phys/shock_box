function [mu,beta,kappa] = artificial(divu,S)
global p rho perx pery del_A del_B gamma

% Get the artificial shear viscosity
c_mu = .02;
tmpx= four_dx(S,1, perx);
tmpy = four_dx(S,2,pery);
tmp = tmpx.*del_A.^2 + tmpy.*del_B.^2;
tmp = filters(abs(rho.*tmp),'G');
mu = c_mu*tmp;

% Get the artificial bulk viscosity
c_beta = 1;
tmpx = four_dx(divu,1, perx);
tmpy = four_dx(divu,2,pery);
tmp = tmpx.*del_A.^2 + tmpy.*del_B.^2;
tmp = filters(abs(rho.*tmp),'G');
beta = c_beta*tmp;

% Get the artificial thermal conductivity
c_kappa = .001;
T = p ./ rho;
c = sqrt(T*gamma);
tmpx = four_dx(T,1,perx);
tmpy = four_dx(T,2,pery);
tmp = tmpx.*del_A + tmpy.*del_B;
tmp = tmp .* rho .* c.^3 ./ T.^2; 
tmp = filters(abs(tmp),'G');
kappa = c_kappa*tmp;


end
