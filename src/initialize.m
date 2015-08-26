function initialize(case_num)
global rho u v p e g x_c y_c U nA nB probname gamma prob_const

%  Initialize prim-variables
rho = ones(nA,nB);
u = zeros(nA,nB);
v = zeros(nA,nB);
p = ones(nA,nB);
e = ones(nA,nB);

switch probname
        
    case{'advect1d'}
        h0 = 1.0;
        amp = 5;
        width = pi/16;
        cenx = max( max( x_c) )/2;
        rad = sqrt( (x_c-cenx).^2 );
        rho = amp*h0*exp( -rad.^2/(width) ) + h0;
        p = rho;
        
    case{'advect2d_cart'}
        h0 = 1.0;
        amp = 5;
        width = pi/16;
        cenx = max( max( x_c) )/2;
        ceny = max( max( y_c) )/2;
        rad = sqrt( (x_c-cenx).^2 + (y_c -ceny).^2);
        rho = amp*h0*exp( -rad.^2/(width) ) + h0;
        dtheta = 2*pi/(nA);
        p = rho;
    
    case{'advect2d_wavy'}
        h0 = 1.0;
        amp = 5;
        width = pi/16;
        cenx = max( max( x_c) )/2;
        ceny = max( max( y_c) )/2;
        rad = sqrt( (x_c-cenx).^2 + (y_c -ceny).^2);
        rho = amp*h0*exp( -rad.^2/(width) ) + h0;
        dtheta = 2*pi/(nA);
        p = rho;
        
    case{'gauss2d_cart'}
        h0 = 1.0;
        amp = 5;
        width = pi/32;
        cenx = max( max( x_c) )/2;
        ceny = max( max( y_c) )/2;
        rad = sqrt( (x_c-cenx).^2 + (y_c -ceny).^2);
        rho = amp*h0*exp( -rad.^2/(width) ) + h0;
        dtheta = 2*pi/(nA);
        %rho = 1-.5 * ( tanh( (x_c - cenx)/ dtheta  )  + 1 );
        %rho = h0 + amp*h0*rho;
        p = rho;
    
    case{'gauss2d_wavy'}
        h0 = 1.0;
        amp = 5;
        width = pi/32;
        cenx = max( max( x_c) )/2;
        ceny = max( max( y_c) )/2;
        rad = sqrt( (x_c-cenx).^2 + (y_c -ceny).^2);
        rho = amp*h0*exp( -rad.^2/(width) ) + h0;
        dtheta = 2*pi/(nA);
        p = rho;
        
    case{'cylinder'}
        h0 = 1.0;
        rho(:,:) = h0;
        p = rho;
        u(:,:) = 2.0 * sqrt(gamma);
                
    case{'bullet'}
        h0 = 1.0;
        rho(:,:) = h0;
        p = rho;
        u(:,:) = prob_const * sqrt(gamma);
        
    case{'shocktube'}
        h0 = 1.0;
        amp = 5;
        cenx = max( max( x_c) )/2;;
        dtheta = 2*pi/(nA);
        rho = 1-.5 * ( tanh( (x_c - cenx)/ dtheta  )  + 1 );
        rho = h0 + amp*h0*rho;
        p = rho;
        
        
    case{'ramp'}
        h0 = 1.0;
        rho(:,:) = h0;
        p = rho;
        u(:,:) = 2.0 * sqrt(gamma);
        
            
    case{'nozzle'}
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
        %get_grid(nA,nB,1,1);
        
        
end

e = p./(gamma-1) + .5 * rho.*(u.*u + v.*v);
U(:,:,1) = rho;
U(:,:,2) = rho.*u;
U(:,:,3) = rho.*v;
U(:,:,4) = e;


end