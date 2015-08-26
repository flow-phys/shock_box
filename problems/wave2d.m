function wave2d(method);
global rho u v p e U probname nA nB r perx pery gamma x_c y_c prob_const grid iter dt tstep

switch method
    
    case    {'setup'}    
        grid = 'cart';
        nA = 64;               %  Number of grid cells in Xchi direction 
        nB = 64;                 %  Number of grid cells in Eta direction
        perx = 1;
        pery = 1;
        iter = 1000;
        dt = 1.0e-2;
        gamma = 1.4;
        if (grid=='cart')
            get_grid(nA,nB,2*pi,2*pi);    % Get grid,metrics, and plot it
        elseif (grid=='wavy')
            amp = .1;
            wave = 1;
            get_grid(nA,nB,amp,wave);
        end
        tstep = 'rk4';
        
    case    {'init'}
        h0 = 1.0;
        amp = 5;
        width = pi/32;
        cenx = max( max( x_c) )/2;
        ceny = max( max( y_c) )/2;
        rad = sqrt( (x_c-cenx).^2 + (y_c -ceny).^2);
        rho = amp*h0*exp( -rad.^2/(width) ) + h0;
        p = rho;
        
        %  Put back into conserved form
        e = p/(gamma-1) + .5 * rho.*(u.*u + v.*v);      
        U(:,:,1) = rho;
        U(:,:,2) = rho.*u;
        U(:,:,3) = rho.*v;
        U(:,:,4) = e; 
        
    case    {'bound'}
        
        %  Periodic          
        
    case    {'viz'}
        plotter(1) = 4;
        viz(plotter);

end