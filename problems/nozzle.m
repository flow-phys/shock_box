function nozzle(method);
global rho u v p e U probname nA nB r perx pery gamma x_c y_c prob_const grid iter dt tstep

switch method        
    
    case    {'setup'}
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
        get_grid(nA,nB,gridfile,2);    % Get grid,metrics, and plot itit  
        
    case    {'init'}
        NOZ = load('init.tec');
        count = 1;
        for j=1:nB
            for i=1:nA
                %x_c(i,j) = NOZ(count,1);
                %y_c(i,j) = NOZ(count,2);
                p(i,j) = NOZ(count,3);
                rho(i,j) = NOZ(count+1,1);
                u(i,j) = NOZ(count+1,2);
                v(i,j) = NOZ(count+1,3);  
                count = count + 2;
            end 
        end
        clear NOZ;
        
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
        u(:,1) = u(:,2);
        u(nA,:) = u(nA-1,:);
        u(:,nB) = 0;
        v(:,1) = 0;
        v(1,:) = 0;
        v(nA,:) = v(nA-1,:);
        v(:,nB) = 0;
        rho(:,1) = rho(:,2); 
        rho(1,:) =rho(2,:);
        rho(nA,:) = rho(nA-1,:);
        rho(:,nB) = rho(:,nB-1);
        p(:,1) = p(:,2); 
        p(1,:) =p(2,:);
        p(nA,:) = p(nA-1,:);
        p(:,nB) = p(:,nB-1); 
        
        %  Put back into conserved form
        e = p/(gamma-1) + .5 * rho.*(u.*u + v.*v);      
        U(:,:,1) = rho;
        U(:,:,2) = rho.*u;
        U(:,:,3) = rho.*v;
        U(:,:,4) = e;
                

end