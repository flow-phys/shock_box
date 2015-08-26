function blunt(method);
global rho u v p e U probname nA nB r perx pery gamma x_c y_c prob_const grid iter dt tstep restart cfl
global step

p0 = 1;
rho0 = 1;
mach = 2.2;
r1 = 1;
r2 = 4;


switch method        
    
    case    {'setup'}
        grid = 'prob';
        nA = 64;               %  Number of grid cells in Xchi direction 
        nB = 32;                 %  Number of grid cells in Eta direction
        iter = 600;
        dt = 0.5e-2;
        cfl = 0.6;
        restart = 0;
        tstep = 'rk4';
        R = 4*pi;
        r = pi/2;
        perx = 0;
        pery = 0;
        gamma = 1.4;
        get_grid(nA,nB,R,r);    % Get grid,metrics, and plot it  
       
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
        %plotter(2) = 6;
        viz(plotter);
        
    case    {'bound'}
        
        u(1,:) = u(2,:);
        v(1,:) = 0;
        p(1,:) = p(2,:);
        %bound('An_damp',
        %bound('B1_extrap');
        arg = [ 2 , p0, rho0, mach];
        bound('Bn_damp', arg);
        bound('B1_slip',0);

        
        
        %  Put back into conserved form
        e = p/(gamma-1) + .5 * rho.*(u.*u + v.*v);      
        U(:,:,1) = rho;
        U(:,:,2) = rho.*u;
        U(:,:,3) = rho.*v;
        U(:,:,4) = e;
        
        if (step < 60)
            %  Filter for stability
            for i=1:4
                U(:,:,i) = filters( U(:,:,i) , 'G');
            end
        end    
        
    case    {'grid'}
        dtheta = pi/2  / (nA-1);
        dR = (r2-r1) / (nB-1); 

        for i=1:nA
            for j=1:nB
                R = r1 + (j-1)*dR;
                tht = pi - dtheta*(i-1);
                x_c(i,j) = R * cos(tht);
                y_c(i,j) = R * sin(tht);
            end
        end

end