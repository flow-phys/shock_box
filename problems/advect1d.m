function advect1d(method);
global rho u v p e U probname nA nB r perx pery gamma x_c y_c prob_const grid iter dt tstep
global time rho_i

switch method        
    
    case    {'setup'}
        grid = 'cart';
        nA = 64;               %  Number of grid cells in Xchi direction 
        nB = 1;                 %  Number of grid cells in Eta direction
        Lx = 2*pi;
        Ly = 0;
        speed = 1;  
        prob_const = speed;
        perx = 1;
        pery = 1;
        iter = 1000;
        dt = 1.0e-3;
        gamma = 1.4;
        get_grid(nA,nB,Lx,Ly);    % Get grid,metrics, and plot it
        tstep = 'rk4';   
        
    case    {'init'}
        h0 = 1.0;
        amp = 5;
        width = pi/16;
        cenx = max( max( x_c) )/2;
        rad = sqrt( (x_c-cenx).^2 );
        rho = amp*h0*exp( -rad.^2/(width) ) + h0;
        p = rho;
        rho_i = rho;
        
        %  Put back into conserved form
        e = p/(gamma-1) + .5 * rho.*(u.*u + v.*v);      
        U(:,:,1) = rho;
        U(:,:,2) = rho.*u;
        U(:,:,3) = rho.*v;
        U(:,:,4) = e; 
        
    case    {'viz'}
        plotter(1) = 1;
        viz(plotter);
        %  Plot the analytical solution
        h0 = 1.0;
        amp = 5;
        width = pi/16;
        cenx = max( max( x_c) )/2;
        cenx = cenx + prob_const*time;
        if (cenx > 3*pi)
            cenx = cenx - ( 2*pi - 2*pi /(nA-1))  ; %(nA/(nA-1));
        end
        rad = sqrt( (x_c-cenx).^2 );
        rho_i = amp*h0*exp( -rad.^2/(width) ) + h0;
        cenx = cenx -  ( 2*pi - 2*pi /(nA-1))   ;
        rad = sqrt( (x_c-cenx).^2 );
        rho_i = rho_i + amp*h0*exp( -rad.^2/(width) );
        hold on;plot(x_c(:,:), rho_i(:,:) ,'ro');drawnow;hold off;
        
    case    {'bound'}
        u = prob_const;
        v = 0;
        p = 1;  
        
        %  Put back into conserved form
        e = p/(gamma-1) + .5 * rho.*(u.*u + v.*v);      
        U(:,:,1) = rho;
        U(:,:,2) = rho.*u;
        U(:,:,3) = rho.*v;
        U(:,:,4) = e;
                

end