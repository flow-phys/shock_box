%% get_grid.m- Obtain the physical grid and calculate the metric terms
%% needed for the general coordinate frame transformation
%% Mesh type variable 'grid' is set by prob('setup') or can be explicitly
%% given using a prob('grid') subroutine and setting grid='prob'

function get_grid(nA,nB,p1,p2)
global x_c y_c dr theta grid perx pery gperx gpery

%% By default, grid is assumed to be non-periodic
gperx = 0;
gpery = 0;

switch grid
    
    %% Cylindrical Mesh
    case{'cyl'}
        R = p1;
        r = p2;
        dtheta = 2*pi/(nA);
        dr = (R-r)/(nB-1);
        for i=1:nA
            for j=1:nB
                rtmp = r+dr*(j-1);
                theta(i) = -dtheta*(i-1) ;
                x_c(i,j) = rtmp*cos(theta(i));
                y_c(i,j) = rtmp*sin(theta(i));
            end
        end
        gperx = 1;  % Grid is truly periodic in x (A) direction
        gpery = 0;  % Grid is NOT periodic in y (B) direction
        
    %% Cartesian Mesh
    case{'cart'}
        dx = p1/(nA-1);
        dy = p2/(nB-1);
        for i=1:nA
            for j=1:nB
                x_c(i,j) = (i-1)*dx;
                y_c(i,j) = (j-1)*dy;
            end
        end
        if (nA==1); x_c=p1; end;
        if (nB==1); y_c=p2; end;
    
    %  Read in a specific mesh     
    case{'read'}
        gridfile = p1;
        griddata = load(gridfile);
        count = 1;cp = p2;
        for j=1:nB
            for i=1:nA
                x_c(i,j) = griddata(count,1);
                y_c(i,j) = griddata(count,2);
                count = count + p2;
            end
        end
        clear griddata;
    
    %  Wavy sine wave mesh    
    case{'wavy'}
        dx = 2*pi / (nA - 1);
        dy = dx;
        ampx = p1*(nA*dx);
        ampy = p1*(nB*dy);

        for i = 1 : nA
            for j = 1 : nB
                sinx(i) = ampx*sin(p2*(i-1)*2*pi/nA);
                siny(j) = ampy*sin(p2*(j-1)*2*pi/nB);
                x_c(i,j) = dx*(i-1) + siny(j);
                y_c(i,j) = dy*(j-1) - sinx(i);
            end
        end

        
    case{'prob'}
        %% You can set x_c, y_c gperx, gpery here.
        prob('grid');
    
end

%% Switch to grid periodicity, setup matrices and compute the metrics
tperx = perx;
tpery = pery;
perx = gperx;pery = gpery;
setup_matrices();

get_metrics();

%% Switch back to the problem periodicity and setup matrices.
perx = tperx;
pery = tpery;
setup_matrices();


end

function plot_grid()
global x_c y_c dxdA
figure(3);clf;
hold all;
plot(x_c,y_c,'k')
for i = 1 : size(x_c,1)
    plot(x_c(i,:),y_c(i,:),'k')
end

end

function get_metrics()
global x_c y_c dAdx dAdy dBdx dBdy Dxy r grid perx pery del_A del_B nA nB dA dB

per = 0;
switch grid
    case{'cyl'}
        per = 1;
    case{'cart'}
        per = 0;
end

dxdA = deriv( x_c, 1, perx*per );
dxdB = deriv( x_c, 2, pery*per );
dydA = deriv( y_c, 1, perx*per );
dydB = deriv( y_c, 2, pery*per );

if (nA == 1); dxdA=1; dydA=1; end;
if (nB == 1); dxdB=1; dydB=1; end;

Dxy = dxdA.*dydB - dxdB.*dydA;


dBdy = dxdA./Dxy;
dBdx = -dydA./Dxy;
dAdx = dydB./Dxy;
dAdy = -dxdB./Dxy;

dA = 1;
dB = 1;
del_A = sqrt( (dxdA*dA).^2 + (dydA*dA).^2  );
del_B = sqrt( (dxdB*dB).^2 + (dydB*dB).^2  );


end


