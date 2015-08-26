function shear(method);
global rho u v p e U probname nA nB r perx pery gamma x_c y_c prob_const grid iter dt tstep restart

switch method        
    
    case    {'setup'}
        grid = 'cart';
        nA = 128;               %  Number of grid cells in Xchi direction 
        nB = 64;                 %  Number of grid cells in Eta direction
        iter = 2000;
        dt = 0.1e-2;
        restart = 1;
        tstep = 'euler';
        perx = 1;
        pery = 0;
        gamma = 1.4;
        Lx = 2;
        Ly = 1;
        get_grid(nA,nB,Lx,Ly);    % Get grid,metrics, and plot it  
        
    case    {'init'}
        rho1=1;
        rho2=4;
        p(:,:) = 1;
        mid = max(max(y_c))/2;
        ceny = mid +  sin(x_c/2*2*pi)*mid/nB*3;%   ( rand(size(y_c)) - .5)*(mid/nB)*1;
        tmp = (1 + tanh( (y_c-ceny)/ (4*mid/nB) )) / 2;
        u = tmp-.5;
        u = u * sqrt(gamma);
        
        rho = rho1 + (rho2 - rho1)*tmp;
        
        %  Put back into conserved form
        e = p/(gamma-1) + .5 * rho.*(u.*u + v.*v);      
        U(:,:,1) = rho;
        U(:,:,2) = rho.*u;
        U(:,:,3) = rho.*v;
        U(:,:,4) = e; 
        
        %  Filter for stability
        for i=1:4
            U(:,:,i) = filters( U(:,:,i) , 'C');
        end

        
    case    {'viz'}
        plotter(1) = 5;
        viz(plotter);
        
    case    {'bound'}
        u(:,1) = u(:,2);
        u(:,nB) = u(:,nB-1);
        v(:,1) = 0;
        v(:,nB) = 0;
        rho(:,1) = rho(:,2); 
        p(:,1) = p(:,2);    
        rho(:,nB) = rho(:,nB-1);
        p(:,nB) =p(:,nB-1);
        
        %  Put back into conserved form
        e = p/(gamma-1) + .5 * rho.*(u.*u + v.*v);      
        U(:,:,1) = rho;
        U(:,:,2) = rho.*u;
        U(:,:,3) = rho.*v;
        U(:,:,4) = e;
                

end