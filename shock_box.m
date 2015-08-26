function shock_box
clc;
clear all;
global nA nB g perx pery probname iter dt cfl U restart A1_U time fil_alpha step cfil


% Add path for problems directory
slash = '/';
if ( findstr( getenv('OS'), 'Windows')  > 0 )
    slash = '\';
end
folder = [ pwd, slash, 'problems'];
path(path,folder);
folder = [ pwd, slash, 'src'];
path(path,folder);

%  Problem parameters
%probname = 'airfoil';
probname = 'cylinder2d';
%probname = 'blunt';
%probname = 'shocktube';
%probname = 'ramp';
%probname = 'nozzle';
%probname = 'wave2d';
%probname = 'advect1d';
probname = 'advect2d';
%probname = 'burgers';
%probname = 'bullet';
%probname = 'shear';
%probname = 'double_mach';

cfl = 1.0;
viz_freq = 20;
mov = 0;
restart = 0;
fil_alpha = .495;
cfil = 'C8';

%  Initialize variables
prob('setup');
allocate();

if (restart == 1)
    restart_run(); 
    disp('Reading restart')
else
    prob('init');
end

count = 1;
time = 0;
prob('viz');
for i=1: iter
    step = i;
    %  Echo time step and take it
    rk4_step();
    disp(['Time step: ',int2str(int16(i))  , '   dt=', num2str(dt) ]);
    %disp(i)
    if ( nandetect() == true )
        disp('Nans detected: shut it down');
        break;
    end
    if (mod(i,viz_freq) == 0)
        count = count +1;
        prob('viz');
        if (mov ==1)
            name = 1000 + i;
            sname = int2str(name);
            saveas(gcf, sname,'tiff')
        end

    end
    
end
prob('viz');
resname = [probname,'.mat'];
save(resname, 'U');

end 

function restart_run()
global  U restart rho u v p e gamma nA nB probname
    
    resname = [probname,'.mat'];
    load (resname, 'U');
    rho = U(:,:,1);
    u = U(:,:,2) ./ rho;
    v = U(:,:,3) ./ rho;
    e = U(:,:,4);
    p = (gamma-1)*(e-.5*rho.*(u.*u+v.*v));

end

function allocate()
global rho u v p e nA nB beta
global U

% Initialize the primative variables
rho = ones(nA,nB);  % Density
u = zeros(nA,nB);    % X-velocity
v = zeros(nA,nB);    % Y-velocity
p = ones(nA,nB);     % Pressure
e = ones(nA,nB);     % Total Energy

U = ones(nA,nB,4);  % Conserved variables

beta = ones(nA,nB);% Artificial Bulk Viscosity

end

function nd = nandetect()
global rho

nd = max(max(isnan(rho)));

end