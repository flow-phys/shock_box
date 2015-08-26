function f = filters(f,type)
global perx pery nA nB

%% X-direction
if (nA ~= 1) 
    switch type
        case{'E5'}
            f = filt_e5_1(f,perx);
        case{'C8'}
            f = filt_c8_1(f,perx);
        case{'G'}
            f = filt_g_1(f,perx);
    end
end

%% Y-direction
if (nB ~= 1)
    switch type
        case{'E5'}
            f = filt_e5_2(f,perx);
        case{'C8'}
            f = filt_c8_2(f,perx);
        case{'G'}
            f = filt_g_2(f,perx);
    end
end


end


function f = filter_c(f,dir,per)
[nx,ny] = size(f);

% Get x-dir
if (dir==1)
n = ny;
for i=1:n
    f(:,i) = filt_c( f(:,i), per );
end
end

if (dir==2)
n = nx;
% Get y-dir
for i=1:n
    f(i,:) = filt_c( f(i,:), per );
end
end

end

function F = filt_c(F,per)
n = max (size(F) );

a1 = 0.6875;                                         
b2 = 0.234375;                                              
c2 = -0.09375;                                               
d2 = 0.015625;

if (per == 1)
    I = [1:n];
    Ip = I + 1;Ip(n)=1;
    Im = I -1; Im(1) = n;
    Ipp = Ip +1; Ipp(n-1)=1;
    Imm = Im -1; Imm(2)=n;
    Ippp = Ipp + 1; Ippp(n-2) = 1;
    Immm = Imm-1;Immm(3) = n;
    
    F(I) = a1* F(I) + b2*( F(Ip) + F(Im) )  + c2 * (F(Ipp) + F(Imm) ) + d2*(F(Ippp) + F(Immm));
    
else
    I = [4:n-3];
    Ip = I + 1;
    Im = I -1;
    Ipp = Ip +1;
    Imm = Im -1; 
    Ippp = Ipp + 1;
    Immm = Imm-1;
    
    F(I) = a1* F(I) + b2*( F(Ip) + F(Im) )  + c2 * (F(Ipp) + F(Imm) ) + d2*(F(Ippp) + F(Immm));
    
    F(2) = .5*F(2) + .25*(F(1)+F(3));
    F(3) = .8125*F(3) + .125*(F(2)+F(4)) + -.03125*(F(1)+F(5));
    F(n-2) = .8125*F(n-2) + .125*(F(n-1)+F(n-3)) + -.03125*(F(n)+F(n-4));
    F(n-1) = .5*F(n-1) + .25*(F(n)+F(n-2));

end
end

function F = filt_e5_1(F,per)
n = max (size(F,1) );

a1 = 0.6875;                                         
b2 = 0.234375;                                              
c2 = -0.09375;                                               
d2 = 0.015625;

if (per == 1)
    I = [1:n];
    Ip = I + 1;Ip(n)=1;
    Im = I -1; Im(1) = n;
    Ipp = Ip +1; Ipp(n-1)=1;
    Imm = Im -1; Imm(2)=n;
    Ippp = Ipp + 1; Ippp(n-2) = 1;
    Immm = Imm-1;Immm(3) = n;
    
    F(I,:) = a1* F(I,:) + b2*( F(Ip,:) + F(Im,:) )  + c2 * (F(Ipp,:) + F(Imm,:) ) + d2*(F(Ippp,:) + F(Immm,:));
    
else
    I = [4:n-3];
    Ip = I + 1;
    Im = I -1;
    Ipp = Ip +1;
    Imm = Im -1; 
    Ippp = Ipp + 1;
    Immm = Imm-1;
    
    F(I,:) = a1* F(I,:) + b2*( F(Ip,:) + F(Im,:) )  + c2 * (F(Ipp,:) + F(Imm,:) ) + d2*(F(Ippp,:) + F(Immm,:));
    
    F(2,:) = .5*F(2,:) + .25*(F(1,:)+F(3,:));
    F(3,:) = .8125*F(3,:) + .125*(F(2,:)+F(4,:)) + -.03125*(F(1,:)+F(5,:));
    F(n-2,:) = .8125*F(n-2,:) + .125*(F(n-1,:)+F(n-3,:)) + -.03125*(F(n,:)+F(n-4,:));
    F(n-1,:) = .5*F(n-1,:) + .25*(F(n,:)+F(n-2,:));

end
end

function F = filt_e5_2(F,per)
n = max (size(F,2) );

a1 = 0.6875;                                         
b2 = 0.234375;                                              
c2 = -0.09375;                                               
d2 = 0.015625;


if (per == 1)
    I = [1:n];
    Ip = I + 1;Ip(n)=1;
    Im = I -1; Im(1) = n;
    Ipp = Ip +1; Ipp(n-1)=1;
    Imm = Im -1; Imm(2)=n;
    Ippp = Ipp + 1; Ippp(n-2) = 1;
    Immm = Imm-1;Immm(3) = n;
    
    F(:,I) = a1* F(:,I) + b2*( F(:,Ip) + F(:,Im) )  + c2 * (F(:,Ipp) + F(:,Imm) ) + d2*(F(:,Ippp) + F(:,Immm));
     
else
    I = [4:n-3];
    Ip = I + 1;
    Im = I -1;
    Ipp = Ip +1;
    Imm = Im -1; 
    Ippp = Ipp + 1;
    Immm = Imm-1;
    
    F(:,I) = a1* F(:,I) + b2*( F(:,Ip) + F(:,Im) )  + c2 * (F(:,Ipp) + F(:,Imm) ) + d2*(F(:,Ippp) + F(:,Immm));
    
    F(:,2) = .5*F(:,2) + .25*(F(:,1)+F(:,3));
    F(:,3) = .8125*F(:,3) + .125*(F(:,2)+F(:,4)) + -.03125*(F(:,1)+F(:,5));
    F(:,n-2) = .8125*F(:,n-2) + .125*(F(:,n-1)+F(:,n-3)) + -.03125*(F(:,n)+F(:,n-4));
    F(:,n-1) = .5*F(:,n-1) + .25*(F(:,n)+F(:,n-2));

end
end


%% 8th-order low pass filter
function dF = filt_c8_1(F,per)
global B1_L B1_U fil_alpha

n = max (size(F,1) ); 

a0 = (93+70*fil_alpha)/128;
a1 = (7+18*fil_alpha)/16;
a2 = (-7+14*fil_alpha)/32;
a3 = (1-2*fil_alpha)/16;
a4 = (-1+2*fil_alpha)/128;

rhs = F;
if (per == 1)
    I = [1:n];
    Ip = I + 1;Ip(n)=1;
    Im = I -1; Im(1) = n;
    Ipp = Ip +1; Ipp(n-1)=1;
    Imm = Im -1; Imm(2)=n;
    Ippp = Ipp +1; Ippp(n-2)=1;
    Immm = Imm -1; Immm(3)=n;
    Ipppp = Ippp +1; Ipppp(n-3)=1;
    Immmm = Immm -1; Immmm(4)=n;

    rhs(I,:) = a0*F(I,:) + a1/2*( F(Ip,:) + F(Im,:) ) + a2/2* (F(Ipp,:) + F(Imm,:)  ) +...
         a3/2* (F(Ippp,:) + F(Immm,:)  ) + a4/2* (F(Ipppp,:) + F(Immmm,:)  );
    
else
    I = [5:n-4];
    Ip = I + 1;
    Im = I -1;
    Ipp = Ip +1;
    Imm = Im -1;
    Ippp = Ipp + 1;
    Immm = Imm - 1;
    Ipppp = Ippp + 1;
    Immmm = Immm - 1;
    
    rhs(I,:) = a0*F(I,:) + a1/2*( F(Ip,:) + F(Im,:) ) + a2/2* (F(Ipp,:) + F(Imm,:)  ) +...
         a3/2* (F(Ippp,:) + F(Immm,:)  ) + a4/2* (F(Ipppp,:) + F(Immmm,:) );
    
    rhs(1,:) = F(1,:);
    rhs(2,:) = .5*F(2,:) + .25*(F(1,:)+F(3,:));
    rhs(3,:) = .8125*F(3,:) + .125*(F(2,:)+F(4,:)) + -.03125*(F(1,:)+F(5,:));
    rhs(4,:) = .8125*F(4,:) + .125*(F(3,:)+F(5,:)) + -.03125*(F(2,:)+F(6,:));
    
    rhs(n-3,:) = .8125*F(n-3,:) + .125*(F(n-2,:)+F(n-4,:)) + -.03125*(F(n-1,:)+F(n-3,:));
    rhs(n-2,:) = .8125*F(n-2,:) + .125*(F(n-1,:)+F(n-3,:)) + -.03125*(F(n,:)+F(n-4,:));
    rhs(n-1,:) = .5*F(n-1,:) + .25*(F(n,:)+F(n-2,:)); 
    rhs(n,:) = F(n,:);
     
end

    dF = B1_U\(B1_L\rhs);

    
end

function dF = filt_c8_2(F,per)
global B2_L B2_U fil_alpha

n = max (size(F,2) ); 

a0 = (93+70*fil_alpha)/128;
a1 = (7+18*fil_alpha)/16;
a2 = (-7+14*fil_alpha)/32;
a3 = (1-2*fil_alpha)/16;
a4 = (-1+2*fil_alpha)/128;

rhs = F;
if (per == 1)
    I = [1:n];
    Ip = I + 1;Ip(n)=1;
    Im = I -1; Im(1) = n;
    Ipp = Ip +1; Ipp(n-1)=1;
    Imm = Im -1; Imm(2)=n;
    Ippp = Ipp +1; Ippp(n-2)=1;
    Immm = Imm -1; Immm(3)=n;
    Ipppp = Ippp +1; Ipppp(n-3)=1;
    Immmm = Immm -1; Immmm(4)=n;

    rhs(:,I) = a0*F(:,I) + a1/2*( F(:,Ip) + F(:,Im) ) + a2/2* (F(:,Ipp) + F(:,Imm)  ) +...
         a3/2* (F(:,Ippp) + F(:,Immm)  ) + a4/2* (F(:,Ipppp) + F(:,Immmm)  );
    
else
    I = [5:n-4];
    Ip = I + 1;
    Im = I -1;
    Ipp = Ip +1;
    Imm = Im -1;
    Ippp = Ipp + 1;
    Immm = Imm - 1;
    Ipppp = Ippp + 1;
    Immmm = Immm - 1;
    
    rhs(:,I) = a0*F(:,I) + a1/2*( F(:,Ip) + F(:,Im) ) + a2/2* (F(:,Ipp) + F(:,Imm)  ) +...
         a3/2* (F(:,Ippp) + F(:,Immm)  ) + a4/2* (F(:,Ipppp) + F(:,Immmm)  ); 
    
    rhs(:,1) = F(:,1);
    rhs(:,2) = .5*F(:,2) + .25*(F(:,1)+F(:,3));
    rhs(:,3) = .8125*F(:,3) + .125*(F(:,2)+F(:,4)) + -.03125*(F(:,1)+F(:,5));
    rhs(:,4) = .8125*F(:,4) + .125*(F(:,3)+F(:,5)) + -.03125*(F(:,2)+F(:,6));
    
    rhs(:,n-3) = .8125*F(:,n-3) + .125*(F(:,n-2)+F(:,n-4)) + -.03125*(F(:,n-1)+F(:,n-3));
    rhs(:,n-2) = .8125*F(:,n-2) + .125*(F(:,n-1)+F(:,n-3)) + -.03125*(F(:,n)+F(:,n-4));
    rhs(:,n-1) = .5*F(:,n-1) + .25*(F(:,n)+F(:,n-2)); 
    rhs(:,n) = F(:,n);
     
end
    rhs2 = transpose(rhs);
    dF2 = B2_U\(B2_L\rhs2);
    dF = transpose(dF2);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Gaussian filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dF = filt_g_1(F,per)

n = max (size(F,1) ); 

a0 = 3565 / 10368;
a1 = 3091 / 12960;
a2 = 1997 / 25920;
a3 = 149 / 12960;
a4 = 107 / 103680;

rhs = F;
if (per == 1)
    I = [1:n];
    Ip = I + 1;Ip(n)=1;
    Im = I -1; Im(1) = n;
    Ipp = Ip +1; Ipp(n-1)=1;
    Imm = Im -1; Imm(2)=n;
    Ippp = Ipp +1; Ippp(n-2)=1;
    Immm = Imm -1; Immm(3)=n;
    Ipppp = Ippp +1; Ipppp(n-3)=1;
    Immmm = Immm -1; Immmm(4)=n;

    rhs(I,:) = a0*F(I,:) + a1*( F(Ip,:) + F(Im,:) ) + a2* (F(Ipp,:) + F(Imm,:)  ) +...
         a3* (F(Ippp,:) + F(Immm,:)  ) + a4* (F(Ipppp,:) + F(Immmm,:)  );
    
else
    I = [5:n-4];
    Ip = I + 1;
    Im = I -1;
    Ipp = Ip +1;
    Imm = Im -1;
    Ippp = Ipp + 1;
    Immm = Imm - 1;
    Ipppp = Ippp + 1;
    Immmm = Immm - 1;
    
    rhs(I,:) = a0*F(I,:) + a1*( F(Ip,:) + F(Im,:) ) + a2* (F(Ipp,:) + F(Imm,:)  ) +...
         a3* (F(Ippp,:) + F(Immm,:)  ) + a4* (F(Ipppp,:) + F(Immmm,:) );
    
    rhs(1,:) = F(1,:);
    rhs(2,:) = .5*F(2,:) + .25*(F(1,:)+F(3,:));
    rhs(3,:) = .8125*F(3,:) + .125*(F(2,:)+F(4,:)) + -.03125*(F(1,:)+F(5,:));
    rhs(4,:) = .8125*F(4,:) + .125*(F(3,:)+F(5,:)) + -.03125*(F(2,:)+F(6,:));
    
    rhs(n-3,:) = .8125*F(n-3,:) + .125*(F(n-2,:)+F(n-4,:)) + -.03125*(F(n-1,:)+F(n-3,:));
    rhs(n-2,:) = .8125*F(n-2,:) + .125*(F(n-1,:)+F(n-3,:)) + -.03125*(F(n,:)+F(n-4,:));
    rhs(n-1,:) = .5*F(n-1,:) + .25*(F(n,:)+F(n-2,:)); 
    rhs(n,:) = F(n,:);
     
end

    dF = rhs;

    
end

function dF = filt_g_2(F,per)
global B1_L B1_U fil_alpha

n = max (size(F,2) ); 

a0 = 3565 / 10368;
a1 = 3091 / 12960;
a2 = 1997 / 25920;
a3 = 149 / 12960;
a4 = 107 / 103680;

rhs = F;
if (per == 1)
    I = [1:n];
    Ip = I + 1;Ip(n)=1;
    Im = I -1; Im(1) = n;
    Ipp = Ip +1; Ipp(n-1)=1;
    Imm = Im -1; Imm(2)=n;
    Ippp = Ipp +1; Ippp(n-2)=1;
    Immm = Imm -1; Immm(3)=n;
    Ipppp = Ippp +1; Ipppp(n-3)=1;
    Immmm = Immm -1; Immmm(4)=n;

    rhs(:,I) = a0*F(:,I) + a1*( F(:,Ip) + F(:,Im) ) + a2* (F(:,Ipp) + F(:,Imm)  ) +...
         a3* (F(:,Ippp) + F(:,Immm)  ) + a4* (F(:,Ipppp) + F(:,Immmm)  );
    
else
    I = [5:n-4];
    Ip = I + 1;
    Im = I -1;
    Ipp = Ip +1;
    Imm = Im -1;
    Ippp = Ipp + 1;
    Immm = Imm - 1;
    Ipppp = Ippp + 1;
    Immmm = Immm - 1;
    
    rhs(:,I) = a0*F(:,I) + a1*( F(:,Ip) + F(:,Im) ) + a2* (F(:,Ipp) + F(:,Imm)  ) +...
         a3* (F(:,Ippp) + F(:,Immm)  ) + a4* (F(:,Ipppp) + F(:,Immmm)  ); 
    
    rhs(:,1) = F(:,1);
    rhs(:,2) = .5*F(:,2) + .25*(F(:,1)+F(:,3));
    rhs(:,3) = .8125*F(:,3) + .125*(F(:,2)+F(:,4)) + -.03125*(F(:,1)+F(:,5));
    rhs(:,4) = .8125*F(:,4) + .125*(F(:,3)+F(:,5)) + -.03125*(F(:,2)+F(:,6));
    
    rhs(:,n-3) = .8125*F(:,n-3) + .125*(F(:,n-2)+F(:,n-4)) + -.03125*(F(:,n-1)+F(:,n-3));
    rhs(:,n-2) = .8125*F(:,n-2) + .125*(F(:,n-1)+F(:,n-3)) + -.03125*(F(:,n)+F(:,n-4));
    rhs(:,n-1) = .5*F(:,n-1) + .25*(F(:,n)+F(:,n-2)); 
    rhs(:,n) = F(:,n);
     
end
    dF = rhs;
end

