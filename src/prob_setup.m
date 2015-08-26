function prob_setup()
global  nA nB perx pery grid  gamma probname iter dt tstep restart prob_const


%  Problem parameters
switch probname
        
    case{'advect1d'}
        grid = 'cart';
        nA = 128;               %  Number of grid cells in Xchi direction 
        nB = 1;                 %  Number of grid cells in Eta direction
        Lx = 2*pi;
        Ly = 0;
        speed = 1;  
        prob_const = speed;
        perx = 1;
        pery = 1;
        iter = 1000;
        dt = 1.0e-2;
        gamma = 1.4;
        get_grid(nA,nB,Lx,Ly);    % Get grid,metrics, and plot it
        tstep = 'rk4';        
        
    case{'advect2d_cart'}
        grid = 'cart';
        nA = 64;               %  Number of grid cells in Xchi direction 
        nB = 64;                 %  Number of grid cells in Eta direction
        perx = 1;
        pery = 1;
        iter = 1000;
        dt = 1.0e-2;
        speed = 1;  
        prob_const = speed;
        gamma = 1.4;
        get_grid(nA,nB,2*pi,2*pi);    % Get grid,metrics, and plot it
        tstep = 'rk4';
            
    case{'advect2d_wavy'}
        grid = 'wavy';
        nA = 64;               %  Number of grid cells in Xchi direction 
        nB = 64;                 %  Number of grid cells in Eta direction
        perx = 1;
        pery = 1;
        amp = .1;
        wave = 1;
        speed = 1;  
        prob_const = speed;
        iter = 1000;
        dt = 1.0e-2;
        gamma = 1.4;
        get_grid(nA,nB,amp,wave);    % Get grid,metrics, and plot it
        tstep = 'rk4';
        
    case{'gauss2d_cart'}
        grid = 'cart';
        nA = 64;               %  Number of grid cells in Xchi direction 
        nB = 64;                 %  Number of grid cells in Eta direction
        perx = 1;
        pery = 1;
        iter = 1000;
        dt = 1.0e-2;
        gamma = 1.4;
        get_grid(nA,nB,2*pi,2*pi);    % Get grid,metrics, and plot it
        tstep = 'rk4';
            
    case{'gauss2d_wavy'}
        grid = 'wavy';
        nA = 64;               %  Number of grid cells in Xchi direction 
        nB = 64;                 %  Number of grid cells in Eta direction
        perx = 1;
        pery = 1;
        amp = .1;
        wave = 1;
        iter = 300;
        dt = 1.0e-2;
        gamma = 1.4;
        get_grid(nA,nB,amp,wave);    % Get grid,metrics, and plot it
        tstep = 'rk4';
        
    case{'cylinder'}
        grid = 'cyl';
        nA = 128;               %  Number of grid cells in Xchi direction 
        nB = 64;                 %  Number of grid cells in Eta direction
        iter = 450;
        dt = 0.5e-2;
        restart = 0;
        tstep = 'euler';
        R = 4*pi;
        r = pi/2;
        perx = 1;
        pery = 0;
        gamma = 1.4;
        get_grid(nA,nB,R,r);    % Get grid,metrics, and plot it
                    
    case{'bullet'}
        grid = 'cart';
        nA = 256;               %  Number of grid cells in Xchi direction 
        nB = 128;                 %  Number of grid cells in Eta direction
        iter = 1800;
        Mach = 2.9;
        prob_const = Mach;
        dt = 0.1e-2;
        restart = 0;
        tstep = 'rk4';
        Lx = 2*pi;
        Ly = pi;
        perx = 1;
        pery = 0;
        gamma = 1.4;
        get_grid(nA,nB,Lx,Ly);    % Get grid,metrics, and plot it
        
    case{'shocktube'}
        grid = 'cart';
        nA = 256;               %  Number of grid cells in Xchi direction 
        nB = 1;                 %  Number of grid cells in Eta direction
        iter = 500;
        dt = 0.2e-2;
        tstep = 'rk4';
        Lx = 2*pi;
        Ly = 0;
        perx = 0;
        pery = 0;
        gamma = 1.4;
        get_grid(nA,nB,Lx,Ly);    % Get grid,metrics, and plot it
        
    case{'ramp'}
        grid = 'read';
        nA = 80;               %  Number of grid cells in Xchi direction 
        nB = 40;                 %  Number of grid cells in Eta direction
        iter = 600;
        dt = 0.2e-2;
        tstep = 'euler';
        gridfile = 'ramp.xy';
        restart = 0;
        perx = 0;
        pery = 0;
        gamma = 1.4;
        get_grid(nA,nB,gridfile,1);    % Get grid,metrics, and plot it
        
                
    case{'nozzle'}
        grid = 'read';
        nA = 512;               %  Number of grid cells in Xchi direction 
        nB = 64;                 %  Number of grid cells in Eta direction
        iter = 3000;
        dt = 1.0e-2;
        tstep = 'euler';
        gridfile = 'init.tec';
        restart = 0;
        perx = 0;
        pery = 0;
        gamma = 1.4;
        get_grid(nA,nB,gridfile,2);    % Get grid,metrics, and plot it
end



end