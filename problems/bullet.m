function bullet(method);
global rho u v p e U probname nA nB r perx pery gamma x_c y_c prob_const grid iter dt tstep

switch method        
    
    case    {'setup'}
        grid = 'cart';
        nA = 256;               %  Number of grid cells in Xchi direction 
        nB = 128;                 %  Number of grid cells in Eta direction
        iter = 1800;
        Mach = 2.9;
        prob_const = Mach;
        dt = 0.01e-2;
        restart = 0;
        tstep = 'rk4';
        Lx = 2*pi;
        Ly = pi;
        perx = 1;
        pery = 0;
        gamma = 1.4;
        get_grid(nA,nB,Lx,Ly);    % Get grid,metrics, and plot it
        
    case    {'init'}
        h0 = 1.0;
        rho(:,:) = h0;
        p = rho;
        u(:,:) = prob_const * sqrt(gamma);
        
        %  Put back into conserved form
        e = p/(gamma-1) + .5 * rho.*(u.*u + v.*v);      
        U(:,:,1) = rho;
        U(:,:,2) = rho.*u;
        U(:,:,3) = rho.*v;
        U(:,:,4) = e; 
        
    case    {'viz'}
        plotter(1) = 5;
        viz(plotter);
        
    case    {'bound'}
        u(:,nB) = u(:,nB-1);
        v(:,nB) = 0;
        p(:,nB) = p(:,nB-1);
        rho(:,nB) = rho(:,nB-1);
		
        u(1:4,:) = prob_const(1)*sqrt(gamma);
        v(1:4,:) = 0;
        p(1:4,:) = 1;
        rho(1:4,:) = 1;
	
        u(:,1) = u(:,2);
        v(:,1) = 0;
        p(:,1) = p(:,2);
        rho(:,1) = rho(:,2);
	
        %   Blocked
        cenx = max( max( x_c) )/2 - pi/2;
        ceny = max( max( y_c) )/2;
        rad = pi/10;

        radius = sqrt( (x_c - cenx).^2 + (y_c - ceny).^2 );
        tmpr =  1 - (1 + tanh( (radius-rad)/(rad/8) ))/2;
        %wall = tmpr

        u = u.*(1-tmpr);
        v = v.*(1-tmpr);

        u(nA,:) = u(nA-1,:);
        v(nA,:) = v(nA-1,:);
        p(nA,:) = p(nA-1,:);
        rho(nA,:) = rho(nA-1,:);
        
        
        %  Put back into conserved form
        e = p/(gamma-1) + .5 * rho.*(u.*u + v.*v);      
        U(:,:,1) = rho;
        U(:,:,2) = rho.*u;
        U(:,:,3) = rho.*v;
        U(:,:,4) = e;
                

end