function [ddx,ddy] = grad(f,perx,pery)
global dAdx dAdy dBdx dBdy r

ddA = deriv(f,1,perx);
ddB = deriv(f,2,pery);

ddx = ddA.*dAdx + ddB.*dBdx;
ddy = ddA.*dAdy + ddB.*dBdy;

end