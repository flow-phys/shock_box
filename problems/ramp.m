function ramp(method);
global rho u v p e U probname nA nB r perx pery gamma x_c y_c prob_const grid iter dt tstep restart

p0 = 1;
rho0 = 1;
mach = 2;


switch method        
    
    case    {'setup'}
        grid = 'read';
        nA = 80;               %  Number of grid cells in Xchi direction 
        nB = 40;                 %  Number of grid cells in Eta direction
        iter = 100;
        dt = 0.2e-2;
        tstep = 'rk4';
        gridfile = 'ramp.xy';
        restart = 1;
        perx = 0;
        pery = 0;
        gamma = 1.4;
        get_grid(nA,nB,gridfile,1);    % Get grid,metrics, and plot it
        
    case    {'init'}
        rho(:,:) = rho0;
        p = p0;
        u(:,:) = mach * sqrt(gamma*p0/rho0);
        
        
        %  Put back into conserved form
        e = p/(gamma-1) + .5 * rho.*(u.*u + v.*v);      
        U(:,:,1) = rho;
        U(:,:,2) = rho.*u;
        U(:,:,3) = rho.*v;
        U(:,:,4) = e; 
        
    case    {'viz'}
        plotter(1) = 5;
        plotter(2) = 6;
        viz(plotter);
        
    case    {'bound'}
        rho(1,:) = rho0;
        p(1,:) = p0;
        u(1,:) = mach * sqrt(p0*gamma/rho0);
        v(1,:) = 0;
%         Lx = 1;
%         Ly = .2;
%         L = sqrt(Lx^2 + Ly^2);
%         umag =sqrt(  u(:,2).^2 + v(:,2).^2 );
%         
%         utmp = u(1,1) - (u(1,1) - umag*(Lx/L)).* ( tanh(y_c(:,1)*200)  ); 
%         vtmp = v(1,1) - (v(1,1) -  umag*(Ly/L)).* ( tanh(y_c(:,1)*200)  ); 
%         
%         u(:,1) = utmp;
         u(nA,:) = u(nA-1,:);
%         u(:,nB) = u(:,nB-1);
%         v(:,1) = vtmp;
%         v(1,:) = 0;
         v(nA,:) = v(nA-1,:);
%         v(:,nB) = 0;
%         rho(:,1) = rho(:,2); 
         rho(nA,:) = rho(nA-1,:);
%         p(:,1) = p(:,2); 
         p(nA,:) = p(nA-1,:);        

        bound('B1_slip');
        bound('Bn_slip');
        arg = [5,p0,rho0,mach];
        bound('A1_damp',arg);
        bound('An_damp',arg);
        %bound('An_extrap');
            
        %  Put back into conserved form
        e = p/(gamma-1) + .5 * rho.*(u.*u + v.*v);      
        U(:,:,1) = rho;
        U(:,:,2) = rho.*u;
        U(:,:,3) = rho.*v;
        U(:,:,4) = e;
                

end