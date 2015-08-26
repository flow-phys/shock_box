function double_mach(method);
global rho u v p e U probname nA nB r perx pery gamma x_c y_c prob_const grid iter dt tstep restart
global time

switch method        
    
    case    {'setup'}
        grid = 'cart';
        nA = 192;               %  Number of grid cells in Xchi direction 
        nB = 64;                 %  Number of grid cells in Eta direction
        iter = 900;
        dt = 0.01e-2;
        restart = 0;
        tstep = 'euler';
        perx = 0;
        pery = 0;
        gamma = 1.4;
        Lx = 3;
        Ly = 1;
        get_grid(nA,nB,Lx,Ly);    % Get grid,metrics, and plot it  
        %bc_param();
        
    case    {'init'}
        % Mach 10- 60 degrees
        M1 = 10;
        [M2,p2,rho2] = shock_jump(M1,gamma);
        ramp = 45;
        x1 = 1/6;
        rampr = ramp*pi/180;
        rho1=1.4;
        p1 = 1.0;
        rho2=rho1 *  rho2;
        p2 = p1 * p2;
        vps = M1*sqrt(p1*gamma/rho1) - M2*sqrt(p2*gamma/rho2);
        u2 = vps*sin(rampr);
        v2 = -vps*cos(rampr);
        thick = 2/nA * 3;
        
        sqrt (u2^2 + v2^2) / sqrt (p2*gamma/rho2)
        prob_const = [u2,v2,p2,rho2,x1,rampr,thick,M1,M2];
        
        b = x1;  
        m = 1/tan(  rampr );  %run/rise 
        shock = (x_c - m*y_c) - b;

        tmp = (1 + tanh(   shock/ thick )) / 2;
        u = (1-tmp)*u2;
        v = (1-tmp)*v2;
        p = p2 + (p1 - p2)*tmp;
        rho = rho2 + (rho1 - rho2)*tmp;
        
        %  Put back into conserved form
        e = p/(gamma-1) + .5 * rho.*(u.*u + v.*v);      
        U(:,:,1) = rho;
        U(:,:,2) = rho.*u;
        U(:,:,3) = rho.*v;
        U(:,:,4) = e; 
        
        %  Filter for stability
        %for i=1:4
        %    U(:,:,i) = filters( U(:,:,i) , 'C');
        %end

        
    case    {'viz'}
        plotter(1) = 5;
        viz(plotter);
        %figure(2);
        %plot(x_c(:,nB/2),u(:,nB/2)); 
        %hold on;
        %plot(x_c(:,nB/2),v(:,nB/2)); 
        %drawnow;hold off;
        
    case    {'bound'}
        
        u2 = prob_const(1);
        v2 = prob_const(2);
        p2 = prob_const(3);
        rho2 = prob_const(4);
        x1 = prob_const(5);
        rampr = prob_const(6);
        thick = prob_const(7);
        M1 = prob_const(8);
        M2 = prob_const(9);
        
        % Inflow
        u(1,:) = u2;
        v(1,:) = v2;
        rho(1,:) = rho2;
        p(1,:) = p2;
        
        % Bottom reflected wall
        v(:,1) = 0;
        tmp = (1+tanh( (x_c(:,1)-x1)/thick) )/2;
        v(:,1) = v2 + (v(:,1)  - v2 ).*tmp;
        u(:,1) = u2 + (u(:,1)  - u2 ).*tmp;
        p(:,1) = p2 + (p(:,1)  - p2 ).*tmp;
        rho(:,1) = rho2 + (rho(:,1)  - rho2 ).*tmp;
        
        % Top Wall
        sloc = x1 + 1/tan(rampr);
        sloc = sloc + time*M1*sin(rampr)
        %hold on;plot(sloc,1,'ro');drawnow;hold off;
        tmp = (1+tanh( (x_c(:,nB)-sloc)/ thick) )/2;
        v(:,nB) = v2 + ( - v2 ).*tmp;
        u(:,nB) = u2 + (  - u2 ).*tmp;
        p(:,nB) = p2 + (p(:,nB)  - p2 ).*tmp;
        rho(:,nB) = rho2 + (rho(:,nB)  - rho2 ).*tmp;
        
        %  Put back into conserved form
        e = p/(gamma-1) + .5 * rho.*(u.*u + v.*v);      
        U(:,:,1) = rho;
        U(:,:,2) = rho.*u;
        U(:,:,3) = rho.*v;
        U(:,:,4) = e;
                
end
end

function bc_param()
global prob_const gamma nA

        % Mach 10- 60 degrees
        M1 = 10;
        [M2,p2,rho2] = shock_jump(M1,gamma);
        ramp = 60;
        x1 = 1/6;
        rampr = ramp*pi/180;
        rho1=1.4;
        p1 = 1.0;
        rho2=rho1 *  rho2;
        p2 = p1 * p2;
        vps = M1*sqrt(p1*gamma/rho1) - M2*sqrt(p2*gamma/rho2);
        u2 = vps*sin(rampr);
        v2 = -vps*cos(rampr);
        thick = 5/nA * 3;
        
        sqrt (u2^2 + v2^2) / sqrt (p2*gamma/rho2)
        prob_const = [u2,v2,p2,rho2,x1,rampr,thick,M1,M2];
end

function [M2,p2,rho2] = shock_jump(M1,gamma)

    p2 = (2*gamma*M1^2 - (gamma-1)) / (gamma+1);
    rho2 = (gamma+1)*M1^2 / ( (gamma-1)*M1^2 + 2 );
    M2 = sqrt( ((gamma-1)*M1^2 + 2 ) / (2*gamma*M1^2 - (gamma-1)) );

end