function ddf = four_dx(f,dir,per)
[nx,ny] = size(f);

% Get dfdx
if (dir==1)
    if (nx==1)
        ddf = zeros(nx,ny);
    else
        ddf = dfour_1( f,per);
    end
end

if (dir==2)
    if (ny==1)
        ddf = zeros(nx,ny);
    else
        ddf = dfour_2( f, per);
    end
end

end

function ddF = dfour(F,per)
n = max (size(F) );

if (per == 1)
    I = [1:n];
    Ip = I + 1;Ip(n)=1;
    Im = I -1; Im(1) = n;
    Ipp = Ip +1; Ipp(n-1)=1;
    Imm = Im -1; Imm(2)=n;

    ddF(I) = 6*F(I) - 4*( F(Ip) + F(Im) ) + (F(Ipp) + F(Imm)  );
    
else
    I = [3:n-2];
    Ip = I + 1;
    Im = I -1;
    Ipp = Ip +1;
    Imm = Im -1;

    ddF(I) = 6*F(I) - 4*( F(Ip) + F(Im) ) + (F(Ipp) + F(Imm)  );
    
    ddF(1) = ddF(3);
    ddF(2) = ddF(3);
    ddF(n-1) = ddF(n-2);    
    ddF(n) =    ddF(n-2);
    
end


end

function ddF = dfour_1(F,per)
n = max (size(F,1) );

if (per == 1)
    I = [1:n];
    Ip = I + 1;Ip(n)=1;
    Im = I -1; Im(1) = n;
    Ipp = Ip +1; Ipp(n-1)=1;
    Imm = Im -1; Imm(2)=n;

    ddF(I,:) = 6*F(I,:) - 4*( F(Ip,:) + F(Im,:) ) + (F(Ipp,:) + F(Imm,:)  );
    
else
    I = [3:n-2];
    Ip = I + 1;
    Im = I -1;
    Ipp = Ip +1;
    Imm = Im -1;

    ddF(I,:) = 6*F(I,:) - 4*( F(Ip,:) + F(Im,:) ) + (F(Ipp,:) + F(Imm,:)  );
    
    ddF(1,:) = ddF(3,:);
    ddF(2,:) = ddF(3,:);
    ddF(n-1,:) = ddF(n-2,:);    
    ddF(n,:) =    ddF(n-2,:);
    
end


end

function ddF = dfour_2(F,per)
n = max (size(F,2) );

if (per == 1)
    I = [1:n];
    Ip = I + 1;Ip(n)=1;
    Im = I -1; Im(1) = n;
    Ipp = Ip +1; Ipp(n-1)=1;
    Imm = Im -1; Imm(2)=n;

    ddF(:,I) = 6*F(:,I) - 4*( F(:,Ip) + F(:,Im) ) + (F(:,Ipp) + F(:,Imm)  );
    
else
    I = [3:n-2];
    Ip = I + 1;
    Im = I -1;
    Ipp = Ip +1;
    Imm = Im -1;

    ddF(:,I) = 6*F(:,I) - 4*( F(:,Ip) + F(:,Im) ) + (F(:,Ipp) + F(:,Imm)  );
    
    ddF(:,1) = ddF(:,3);
    ddF(:,2) = ddF(:,3);
    ddF(:,n-1) = ddF(:,n-2);    
    ddF(:,n) =    ddF(:,n-2);
    
end


end