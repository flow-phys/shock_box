function cylinder2d(method);
global rho u v p e U probname nA nB r perx pery gamma x_c y_c prob_const grid iter dt tstep restart

switch method        
    
    case    {'setup'}
        grid = 'cyl';
        nA = 64;               %  Number of grid cells in Xchi direction 
        nB = 32;                 %  Number of grid cells in Eta direction
        iter = 300;
        dt = 0.5e-2;
        restart = 0;
        tstep = 'rk4'
        R = 4*pi;
        r = pi/2;
        perx = 1;
        pery = 0;
        gamma = 1.4;
        get_grid(nA,nB,R,r);    % Get grid,metrics, and plot it  
        figure(10);
    case    {'init'}
        h0 = 1.0;
        rho(:,:) = h0;
        p = rho;
        u(:,:) = 2.0 * sqrt(gamma);
        
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
        %u(:,1) = 0;
        %v(:,1) = 0;
        %rho(:,1) = rho(:,2); 
        %p(:,1) = p(:,2);    
        
        bound('B1_extrap');
        bound('B1_slip');
        
        %  Put back into conserved form
        e = p/(gamma-1) + .5 * rho.*(u.*u + v.*v);      
        U(:,:,1) = rho;
        U(:,:,2) = rho.*u;
        U(:,:,3) = rho.*v;
        U(:,:,4) = e;
                

end