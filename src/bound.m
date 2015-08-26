function bound(bound_type,arg)
global rho u v p e U  nA nB perx pery gamma 
global dAdx dAdy dBdx dBdy Dxy x_c y_c

switch bound_type
    
    % Viscous Wall Boundary Condition
    case('A1_wall')
        u(1,:) = 0;
        v(1,:) = 0;
    case('An_wall')
        u(nA,:) = 0;
        v(nA,:) = 0;
    case('B1_wall')
        u(:,1) = 0;
        v(:,1) = 0;
    case('Bn_wall')
        u(:,nB) = 0;
        v(:,nB) = 0;
        
    % Slip Wall (Inviscid) Boundary Condition
    case('A1_slip')
         dydB = dAdx(1,:).*Dxy(1,:);
         dxdB = -dAdy(1,:).*Dxy(1,:);
         t_mag = sqrt( dxdB.^2 + dydB.^2 );     % Magnitude of Tangent vector at wall
         u_mag = sqrt( u(1,:).^2 + v(1,:).^2  );  % Vel. Magnitude
         u(1,:) = u_mag.*dxdB./t_mag;
         v(1,:) = u_mag.*dydB./t_mag;
    case('An_slip')
         dxdA = dBdy(nA,:).*Dxy(nA,:);
         dydA = -dBdx(nA,:).*Dxy(nA,:);
         dydB = dAdx(nA,:).*Dxy(nA,:);
         dxdB = -dAdy(nA,:).*Dxy(nA,:);
         t_mag = sqrt( dxdB.^2 + dydB.^2 );     % Magnitude of Tangent vector at wall
         u_mag = sqrt( u(nA,:).^2 + v(nA,:).^2  );  % Vel. Magnitude
         u(nA,:) = u_mag.*dxdB./t_mag;
         v(nA,:) = u_mag.*dydB./t_mag;
    case('B1_slip')
         dxdA = dBdy(:,1).*Dxy(:,1);
         dydA = -dBdx(:,1).*Dxy(:,1);
         t_mag = sqrt( dxdA.^2 + dydA.^2 );     % Magnitude of Tangent vector at wall
         u_mag = sqrt( u(:,1).^2 + v(:,1).^2  );  % Vel. Magnitude
         sigx = sign(1);
         sigy = sign(dxdA.*dydA);
         u(:,1) = sigx.*abs(u_mag.*dxdA./t_mag);
         v(:,1) = sigy.*abs(u_mag.*dydA./t_mag);
         p(:,1) = p(:,2);
         rho(:,1) = rho(:,2);  
    case('Bn_slip')
         dxdA = dBdy(:,nB).*Dxy(:,nB);
         dydA = -dBdx(:,nB).*Dxy(:,nB);
         t_mag = sqrt( dxdA.^2 + dydA.^2 );     % Magnitude of Tangent vector at wall
         u_mag = sqrt( u(:,nB).^2 + v(:,nB).^2  );  % Vel. Magnitude
         u(:,nB) = u_mag.*dxdA./t_mag;
         v(:,nB) = u_mag.*dydA./t_mag;
         p(:,nB) = p(:,nB-1); 
         rho(:,nB) = rho(:,nB-1);
        
    % Extrapolate the variables
    case('A1_extrap')
         u(1,:) = (54.0*u(2,:)-27.0*u(3,:)+6.0*u(4,:) )/33.0;      
         v(1,:) = (54.0*v(2,:)-27.0*v(3,:)+6.0*v(4,:) )/33.0;  
         p(1,:) = (54.0*p(2,:)-27.0*p(3,:)+6.0*p(4,:) )/33.0;  
         rho(1,:) = (54.0*rho(2,:)-27.0*rho(3,:)+6.0*rho(4,:) )/33.0;  
    case('An_extrap')
         u(nA,:) = (54.0*u(nA-1,:)-27.0*u(nA-2,:)+6.0*u(nA-3,:) )/33.0;      
         v(nA,:) = (54.0*v(nA-1,:)-27.0*v(nA-2,:)+6.0*v(nA-3,:) )/33.0;  
         p(nA,:) = (54.0*p(nA-1,:)-27.0*p(nA-2,:)+6.0*p(nA-3,:) )/33.0;  
         rho(nA,:) = (54.0*rho(nA-1,:)-27.0*rho(nA-2,:)+6.0*rho(nA-3,:) )/33.0;
    case('B1_extrap')
         u(:,1) = (54.0*u(:,2)-27.0*u(:,3)+6.0*u(:,4) )/33.0;      
         v(:,1) = (54.0*v(:,2)-27.0*v(:,3)+6.0*v(:,4) )/33.0;  
         p(:,1) = (54.0*p(:,2)-27.0*p(:,3)+6.0*p(:,4) )/33.0;  
         rho(:,1) = (54.0*rho(:,2)-27.0*rho(:,3)+6.0*rho(:,4) )/33.0;  
    case('Bn_extrap')
         u(:,nB) = (54.0*u(:,nB-1)-27.0*u(:,nB-2)+6.0*u(:,nB-3 ))/33.0;      
         v(:,nB) = (54.0*v(:,nB-1)-27.0*v(:,nB-2)+6.0*v(:,nB-3 ))/33.0;  
         p(:,nB) = (54.0*p(:,nB-1)-27.0*p(:,nB-2)+6.0*p(:,nB-3 ))/33.0;  
         rho(:,nB) = (54.0*rho(:,nB-1)-27.0*rho(:,nB-2)+6.0*rho(:,nB-3 ))/33.0;  
    % Filter the variables near the Boundary
    case('A1_damp')
        npts = arg(1);
        p0 = arg(2);
        rho0 = arg(3);
        mach = arg(4);
        
        %Gradually apply the exit boundary conditions  
        filpt = npts;
        thick = 3;
        for i=1:nA
            dumT(i,1:nB) = (1+tanh( (  i  - filpt) /thick ) ) / 2;
        end
        
        dumT = 1 - dumT;
        % Force back pressure to remain ambient and let density float from NSCBC
        p = p + dumT .*  ( p0 - p ) ;                                           
        u = u + dumT .*  ( mach*sqrt(p0*gamma/rho0) - u ) ; 
        v = v + dumT .*  (  - v ) ; 
        
        % Guassian Filter for last N points in x-direction (A)
        dumF = filters(u,'G');
        u = u + dumT.*(dumF-u);
        dumF = filters(v,'G');
        v = v + dumT.*(dumF-v);
        dumF = filters(rho,'G');
        rho = rho + dumT.*(dumF-rho);
        dumF = filters(p,'G');
        p = p + dumT.*(dumF-p);
      
    case('An_damp')
        npts = arg(1);
        p0 = arg(2);
        rho0 = arg(3);
        mach = arg(4);
        
        %Gradually apply the exit boundary conditions  
        filpt = nA - npts;
        thick = 3;
        for i=1:nA
            dumT(i,1:nB) = (1+tanh( (  i  - filpt) /thick ) ) / 2;
        end
        
        % Force back pressure to remain ambient and let density float from NSCBC
        %p = p + dumT .*  ( p0 - p ) ;                                           
        %u = u + dumT .*  ( mach*sqrt(p0*gamma/rho0) - u ) ; 
        %v = v + dumT .*  (  - v ) ; 
        
        % Guassian Filter for last N points in x-direction (A)
        dumF = filters(u,'G');
        u = u + dumT.*(dumF-u);
        dumF = filters(v,'G');
        v = v + dumT.*(dumF-v);
        dumF = filters(rho,'G');
        rho = rho + dumT.*(dumF-rho);
        dumF = filters(p,'G');
        p = p + dumT.*(dumF-p);
    case('B1_damp')
        
    case('Bn_damp')
        npts = arg(1);
        p0 = arg(2);
        rho0 = arg(3);
        mach = arg(4);
        
        %Gradually apply the exit boundary conditions  
        filpt = nB-npts;
        thick = 3;
        for i=1:nB
            dumT(1:nA,i) = (1+tanh( (  i  - filpt) /thick ) ) / 2;
        end
        
        % Force back pressure to remain ambient and let density float from NSCBC
        p = p + dumT .*  ( p0 - p ) ;                                           
        u = u + dumT .*  ( mach*sqrt(p0*gamma/rho0) - u ) ; 
        v = v + dumT .*  (  - v ) ; 
        
        % Guassian Filter for last N points in x-direction (A)
        dumF = filters(u,'G');
        u = u + dumT.*(dumF-u);
        dumF = filters(v,'G');
        v = v + dumT.*(dumF-v);
        dumF = filters(rho,'G');
        rho = rho + dumT.*(dumF-rho);
        dumF = filters(p,'G');
        p = p + dumT.*(dumF-p);
        
       
        
end

%  Put back into conserved form
e = p/(gamma-1) + .5 * rho.*(u.*u + v.*v);
U(:,:,1) = rho;
U(:,:,2) = rho.*u;
U(:,:,3) = rho.*v;
U(:,:,4) = e;

end