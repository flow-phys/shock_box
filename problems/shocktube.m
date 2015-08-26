function shocktube(method);
global rho u v p e U probname nA nB r perx pery gamma x_c y_c prob_const grid iter dt tstep

switch method        
    
    case    {'setup'}
        grid = 'cart';
        nA = 256;               %  Number of grid cells in Xchi direction 
        nB = 1;                 %  Number of grid cells in Eta direction
        iter = 50;
        dt = 0.1e-2;
        tstep = 'rk4';
        Lx = 2*pi;
        Ly = 2*pi*nB/nA;
        perx = 0;
        pery = 0;
        gamma = 1.4;
        get_grid(nA,nB,Lx,Ly);    % Get grid,metrics, and plot it
        
    case    {'init'}
        h0 = 1.0;
        amp = 5;
        cenx = max( max( x_c) )/2;;
        dtheta = 2*pi/(nA);
        rho = 1-.5 * ( tanh( (x_c - cenx)/ dtheta  )  + 1 );
        rho = h0 + amp*h0*rho;
        p = rho;
        
        %  Put back into conserved form
        e = p/(gamma-1) + .5 * rho.*(u.*u + v.*v);      
        U(:,:,1) = rho;
        U(:,:,2) = rho.*u;
        U(:,:,3) = rho.*v;
        U(:,:,4) = e; 
        
        for i=1:4
            U(:,:,i) = filters(U(:,:,i),'G');
        end
        
    case    {'viz'}
        plotter(1) = 2;
        viz(plotter);
        
    case    {'bound'}
        u(1,:) = 0;
        v(1,:) = 0;
        rho(1,:) = rho(2,:); 
        p(1,:) = p(2,:);
        
        u(end,:) = 0;
        v(end,:) = 0;
        rho(end,:) = rho(end-1,:); 
        p(end,:) = p(end-1,:);
        %  Put back into conserved form
        e = p/(gamma-1) + .5 * rho.*(u.*u + v.*v);      
        U(:,:,1) = rho;
        U(:,:,2) = rho.*u;
        U(:,:,3) = rho.*v;
        U(:,:,4) = e;
                

end