%% prob.m - Simple function for pointing to the foo_prob.f file which
%% defines the various prob methods needed for the solver.  
%% For a n

function prob(method)
global probname

switch probname
    case    {'airfoil'}
        airfoil(method);
    case    {'advect1d'}
        advect1d(method);
    case    {'burgers'}
        burgers(method);
    case    {'blunt'}
        blunt(method);
    case    {'advect2d'}
        advect2d(method);
    case    {'wave2d'}
        wave2d(method);
    case    {'shocktube'}
        shocktube(method);
    case    {'cylinder2d'}
        cylinder2d(method);
    case    {'ramp'}
        ramp(method);
    case    {'nozzle'}
        nozzle(method);
    case    {'bullet'}
        bullet(method);
    case    {'shear'}
        shear(method);
    case    {'double_mach'}
        double_mach(method);
end



end 