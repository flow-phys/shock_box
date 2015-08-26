function rk4_step()
global U F G nA tstep time dt cfl cfil

Ark(1) = 0.0;
Ark(2) = -6234157559845/12983515589748;
Ark(3) = -6194124222391/4410992767914;
Ark(4) = -31623096876824/15682348800105;
Ark(5) = -12251185447671/11596622555746;

Brk(1) = 494393426753/4806282396855;
Brk(2) = 4047970641027/5463924506627;
Brk(3) = 9795748752853/13190207949281;
Brk(4) = 4009051133189/8539092990294;
Brk(5) = 1348533437543/7166442652324;

eta(1) = 494393426753/4806282396855;
eta(2) = 4702696611523/9636871101405;
eta(3) = 3614488396635/5249666457482;
eta(4) = 9766892798963/10823461281321;
eta(5) = 1.0;

%	Initialize some intermediate arrays
tmp1 = zeros(size(U));
tmp2 = zeros(size(U));
PHI = zeros(size(U));

%	Get primative flow variables
get_var();
prob('bound');
dt = CFL(cfl);
%dt = .00001;
time_i = time;
switch tstep
    case{'rk4'}
        for ii=1:5
            %    ii
            get_flux();
            FLUX = diverg(F,G);
            tmp1 =  Ark(ii)*PHI;
            PHI = -dt*FLUX + tmp1;
            tmp2 =  Brk(ii)*PHI;
            U =  U + tmp2; 
            time = time_i + eta(ii)*dt;
            get_var();
            prob('bound');
            
            %  Filter for stability
            for i=1:4
                U(:,:,i) = filters( U(:,:,i) , cfil);
            end
        end

    case{'euler'}
        get_flux();
        FLUX = diverg(F,G);
        U = U - dt*FLUX;
        time = time_i + dt;
        get_var();
        prob('bound');
        %  Filter for stability
        for i=1:4
            U(:,:,i) = filters( U(:,:,i) , cfil);
        end
end



end

function get_var()
global rho u v p e U gamma

    rho = U(:,:,1);
    u = U(:,:,2)./ rho;
    v = U(:,:,3) ./ rho;
    e = U(:,:,4);
    p = (gamma-1)*(e-.5*rho.*(u.*u+v.*v));
    
end


function get_flux()
global rho u v F G g p e
    
    [tauxx,tauxy,tauyy,Qx,Qy] = get_tau();

    F(:,:,1) = rho.*u;
    F(:,:,2) = rho.*u.*u + p - tauxx;
    F(:,:,3) = rho.*u.*v - tauxy;
    F(:,:,4) = (e + p).*u - tauxx.*u - tauxy.*v + Qx;
    
    G(:,:,1) = rho.*v;
    G(:,:,2) = rho.*u.*v - tauxy;
    G(:,:,3) = rho.*v.*v + p - tauyy;
    G(:,:,4) = (e+p).*v - tauxy.*u - tauyy.*v + Qy;
    
end

function [tauxx,tauxy,tauyy,Qx,Qy] = get_tau();
global u v rho p perx pery beta

    mu_f = 0;
    beta_f = 0;
    kappa_f = 0;

    [dudx,dudy] = grad(u,perx,pery);
    [dvdx,dvdy] = grad(v,perx,pery);
    divu = dudx + dvdy;
    S = sqrt( (dudx.^2 + dvdy.^2 + .5*(dudy+dvdx).^2)  );

    [mu_s,beta_s,kappa_s] = artificial(divu,S);
    mu = mu_f + mu_s;
    beta = beta_f  + beta_s;
    kappa = kappa_f  + kappa_s;

    %  Viscous stress
    tauxx = 2*mu.*dudx + ( beta - (2/3)*mu ).*(dudx + dvdy);
    tauxy = mu.*(dudy + dvdx);
    tauyy = 2*mu.*dvdy + ( beta - (2/3)*mu ).*(dudx + dvdy);
    
    %  Heat flux
    T = p./rho;
    [Qx,Qy] = grad(T,perx,pery);
    Qx = -Qx .* kappa;
    Qy = -Qy .* kappa;

end


function RHS = diverg(F,G)
global dAdx dAdy dBdx dBdy Dxy r perx pery nA dr nB
    for i=1:4
        Fh = F(:,:,i) .* dAdx + G(:,:,i) .* dAdy; 
        Gh = F(:,:,i) .* dBdx + G(:,:,i) .* dBdy;  
        dFhdA = deriv( Dxy .* Fh, 1, perx );
        dGhdB = deriv( Dxy .* Gh, 2, pery );
        RHS(:,:,i) = ( dFhdA + dGhdB )./Dxy;
    end
end

function dt = CFL(cfl)
global u v p rho gamma beta
global dAdx dAdy dBdx dBdy Dxy nA nB dA dB del_A del_B 
    
    shock_stab = .25;

    % Transpose jacobian
    dxdA = dBdy.*Dxy;
    dydA = -dBdx.*Dxy;
    dydB = dAdx.*Dxy;
    dxdB = -dAdy.*Dxy;
    u_A = u.*(dxdA*dA./del_A) + v.*(dydA*dA./del_A);
    u_B = u.*(dxdB*dB./del_B) + v.*(dydB*dB./del_B);
    c = sqrt( p*gamma./rho);
    dt = max(max( abs(u_A)./del_A + abs(u_B)./del_B + ...
        c.*sqrt( 1./del_A.^2 + 1./del_B.^2)   ));
    dt_conv = cfl/dt;
    
    del = min(del_A,del_B);
    dt_shock = max(max(beta./(rho.*del.^2)));
    
    if  (dt_shock <= 0) 
        dt_shock = 10.0*dt_conv;
    else
        dt_shock = shock_stab*1.0/dt_shock;
    end
    
    dt = min(dt_conv,dt_shock);
end
