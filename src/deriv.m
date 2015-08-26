function dfdx = deriv(f,dir,per)
[nx,ny] = size(f);

% Get dfdx
if (dir==1)
    if (nx==1)
        dfdx = zeros(nx,ny);
    else
        %dfdx = fourth_1( f, per);
        dfdx = compact6_1(f, per);
        %dfdx = second_1( f, per);
    end
end

if (dir==2)
    if (ny==1)
        dfdx = zeros(nx,ny);
    else
        %dfdx = fourth_2( f, per);
        %dfdx = second_2( f, per);
        dfdx = compact6_2(f, per);
    end
end

end

function dF = fourth(F,per)
n = max (size(F) );

if (per == 1)
    I = [1:n];
    Ip = I + 1;Ip(n)=1;
    Im = I -1; Im(1) = n;
    Ipp = Ip +1; Ipp(n-1)=1;
    Imm = Im -1; Imm(2)=n;

    dF(I) = 2/3*( F(Ip) - F(Im) ) - 1/12* (F(Ipp) - F(Imm)  );
    
else
    I = [3:n-2];
    Ip = I + 1;
    Im = I -1;
    Ipp = Ip +1;
    Imm = Im -1;

    dF(I) = 2/3*( F(Ip) - F(Im) ) - 1/12* (F(Ipp) - F(Imm)  );
    
    dF(1) =-11/6 * F(1) + 3*F(2) - 3/2*F(3) + 1/3*F(4);
    dF(2) = -1/3 *F(1) - 1/2*F(2) + F(3) - 1/6*F(4);
    dF(n-1) = 1/3 *F(n)  + 1/2*F(n-1) - F(n-2)      +  1/6*F(n-3);    
    dF(n) =    11/6 * F(n) - 3*F(n-1)  +  3/2*F(n-2) - 1/3*F(n-3);
    
end

end

function dF = fourth_1(F,per)
n = max (size(F,1) );

if (per == 1)
    I = [1:n];
    Ip = I + 1;Ip(n)=1;
    Im = I -1; Im(1) = n;
    Ipp = Ip +1; Ipp(n-1)=1;
    Imm = Im -1; Imm(2)=n;

    dF(I,:) = 2/3*( F(Ip,:) - F(Im,:) ) - 1/12* (F(Ipp,:) - F(Imm,:)  );
    
else
    I = [3:n-2];
    Ip = I + 1;
    Im = I -1;
    Ipp = Ip +1;
    Imm = Im -1;

    dF(I,:) = 2/3*( F(Ip,:) - F(Im,:) ) - 1/12* (F(Ipp,:) - F(Imm,:)  );
    
    dF(1,:) =-11/6 * F(1,:) + 3*F(2,:) - 3/2*F(3,:) + 1/3*F(4,:);
    dF(2,:) = -1/3 *F(1,:) - 1/2*F(2,:) + F(3,:) - 1/6*F(4,:);
    dF(n-1,:) = 1/3 *F(n,:)  + 1/2*F(n-1,:) - F(n-2,:)      +  1/6*F(n-3,:);    
    dF(n,:) =    11/6 * F(n,:) - 3*F(n-1,:)  +  3/2*F(n-2,:) - 1/3*F(n-3,:);
    
end

end

function dF = fourth_2(F,per)
n = max (size(F,2) );

if (per == 1)
    I = [1:n];
    Ip = I + 1;Ip(n)=1;
    Im = I -1; Im(1) = n;
    Ipp = Ip +1; Ipp(n-1)=1;
    Imm = Im -1; Imm(2)=n;

    dF(:,I) = 2/3*( F(:,Ip) - F(:,Im) ) - 1/12* (F(:,Ipp) - F(:,Imm)  );
    
else
    I = [3:n-2];
    Ip = I + 1;
    Im = I -1;
    Ipp = Ip +1;
    Imm = Im -1;

    dF(:,I) = 2/3*( F(:,Ip) - F(:,Im) ) - 1/12* (F(:,Ipp) - F(:,Imm)  );
    
    dF(:,1) =-11/6 * F(:,1) + 3*F(:,2) - 3/2*F(:,3) + 1/3*F(:,4);
    dF(:,2) = -1/3 *F(:,1) - 1/2*F(:,2) + F(:,3) - 1/6*F(:,4);
    dF(:,n-1) = 1/3 *F(:,n)  + 1/2*F(:,n-1) - F(:,n-2)      +  1/6*F(:,n-3);    
    dF(:,n) =    11/6 * F(:,n) - 3*F(:,n-1)  +  3/2*F(:,n-2) - 1/3*F(:,n-3);
    
end

end

function dF = second_1(F,per)
n = max (size(F,1) );

if (per == 1)
    I = [1:n];
    Ip = I + 1;Ip(n)=1;
    Im = I -1; Im(1) = n;

    dF(I,:) = 1/2*( F(Ip,:) - F(Im,:) )  ;
    
else
    I = [2:n-1];
    Ip = I + 1;
    Im = I -1;

    dF(I,:) = 1/2*( F(Ip,:) - F(Im,:) );
    
    dF(1,:) =-F(1,:) + F(2,:) ; 
    dF(n,:) =    F(n,:) -F(n-1,:);
    
end

end

function dF = second_2(F,per)
n = max (size(F,2) );

if (per == 1)
    I = [1:n];
    Ip = I + 1;Ip(n)=1;
    Im = I -1; Im(1) = n;

    dF(:,I) = 1/2*( F(:,Ip) - F(:,Im) );
    
else
    I = [2:n-1];
    Ip = I + 1;
    Im = I -1;

    dF(:,I) = 1/2*( F(:,Ip) - F(:,Im) ) ;
    
    dF(:,1) =- F(:,1) + F(:,2) ;
    dF(:,n) =   F(:,n) - F(:,n-1) ;
    
end

end


function dF = compact6_1(F,per)
global A1_L A1_U

n = max (size(F,1) ); 

a = 14/9 / 2;
b = 1/9 / 4;
rhs = F;
if (per == 1)
    I = [1:n];
    Ip = I + 1;Ip(n)=1;
    Im = I -1; Im(1) = n;
    Ipp = Ip +1; Ipp(n-1)=1;
    Imm = Im -1; Imm(2)=n;

    rhs(I,:) = a*( F(Ip,:) - F(Im,:) ) + b* (F(Ipp,:) - F(Imm,:)  );
    
else
    I = [3:n-2];
    Ip = I + 1;
    Im = I -1;
    Ipp = Ip +1;
    Imm = Im -1;

    rhs(I,:) = a*( F(Ip,:) - F(Im,:) ) + b* (F(Ipp,:) - F(Imm,:)  );
    
    rhs(1,:) =-11/6 * F(1,:) + 3*F(2,:) - 3/2*F(3,:) + 1/3*F(4,:);
    rhs(2,:) = -1/3 *F(1,:) - 1/2*F(2,:) + F(3,:) - 1/6*F(4,:);
    rhs(n-1,:) = 1/3 *F(n,:)  + 1/2*F(n-1,:) - F(n-2,:)      +  1/6*F(n-3,:);    
    rhs(n,:) =    11/6 * F(n,:) - 3*F(n-1,:)  +  3/2*F(n-2,:) - 1/3*F(n-3,:);
    
end

    dF = A1_U\(A1_L\rhs);

    
end

function dF = compact6_2(F,per)
global A2_L A2_U
n = max (size(F,2) );

a = 14/9 / 2;
b = 1/9 / 4;
rhs = F;

if (per == 1)
    I = [1:n];
    Ip = I + 1;Ip(n)=1;
    Im = I -1; Im(1) = n;
    Ipp = Ip +1; Ipp(n-1)=1;
    Imm = Im -1; Imm(2)=n;

    rhs(:,I) = a*( F(:,Ip) - F(:,Im) )  + b* (F(:,Ipp) - F(:,Imm)  );
    
else
    I = [3:n-2];
    Ip = I + 1;
    Im = I -1;
    Ipp = Ip +1;
    Imm = Im -1;

    rhs(:,I) = a*( F(:,Ip) - F(:,Im) )  + b* (F(:,Ipp) - F(:,Imm)  );
    
    rhs(:,1) =-11/6 * F(:,1) + 3*F(:,2) - 3/2*F(:,3) + 1/3*F(:,4);
    rhs(:,2) = -1/3 *F(:,1) - 1/2*F(:,2) + F(:,3) - 1/6*F(:,4);
    rhs(:,n-1) = 1/3 *F(:,n)  + 1/2*F(:,n-1) - F(:,n-2)      +  1/6*F(:,n-3);    
    rhs(:,n) =    11/6 * F(:,n) - 3*F(:,n-1)  +  3/2*F(:,n-2) - 1/3*F(:,n-3);
    
end
    rhs2 = transpose(rhs);
    dF2 = A2_U\(A2_L\rhs2);
    dF = transpose(dF2);
end


