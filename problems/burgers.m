function burgers(method);
global rho u v p e U probname nA nB r perx pery gamma x_c y_c prob_const grid iter dt tstep
global time rho_i restart

switch method        
    
    case    {'setup'}
        grid = 'cart';
        nA = 256;               %  Number of grid cells in Xchi direction 
        nB = 1;                 %  Number of grid cells in Eta direction
        Lx = 2*pi;
        Ly = 0;
        speed = 1;  
        prob_const = speed;
        perx = 1;
        pery = 1;
        iter = 1000;
        dt = 0.5e-3;
        gamma = 1.4;
        get_grid(nA,nB,Lx,Ly);    % Get grid,metrics, and plot it
        tstep = 'euler';   
        
    case    {'init'}
        h0 = 1.0;
        amp = 5;
        u = amp*sin(2*x_c);
        %u =  - tanh( 10*(x_c - max(max(x_c))/2)) ;
        rho = ones(nA,nB);
        p = -1/2* (u.^2);
        
        %  Put back into conserved form
        e = ones(nA,nB) ;%p/(gamma-1) + .5 * rho.*(u.*u + v.*v);      
        U(:,:,1) = rho;
        U(:,:,2) = rho.*u;
        U(:,:,3) = rho.*v;
        U(:,:,4) = e; 
        
    case    {'viz'}
        plotter(1) = 2;
        viz(plotter);
        
    case    {'bound'}
        u = rho.*u;
        v = zeros(nA,nB);
        p = -rho .*u.*u;
        p = p + 1/2 * u.*u;
        rho = ones(nA,nB);
        
        %  Put back into conserved form
        e = ones(nA,nB);
        U(:,:,1) = rho;
        U(:,:,2) = rho.*u;
        U(:,:,3) = rho.*v;
        U(:,:,4) = e;
                

end