% Sample 2D advection solver on curvilinear coordinates
% This solver is for an example with *radial* coordinates

clear;

% define our curvilinear coordinates (or computational graph) values
mxi = 200; met = 201; % spacing
xi_min = 0.1; xi_max = 1;
et_min = 0; et_max = 2*pi;

% derived parameters
dxi = (xi_max-xi_min)/mxi;
det = (et_max-et_min)/met;

% Grid cell centers
xivec  = linspace(xi_min+0.5*dxi, xi_max-0.5*dxi, mxi );
etvec  = linspace(et_min+0.5*det, et_max-0.5*det, met );

% Grid edges
xievec  = linspace(xi_min, xi_max, mxi+1 );
etevec  = linspace(et_min, et_max, met+1 );

% define our other conditions
Tf = 1.0;                                   % final time
cfl_des = 0.6;                              % desired CFL number
nframes = 10;                               % Number of "frames" to plot

% make our mesh
[XI,ET] = meshgrid(xivec,etvec);
XI = XI'; ET = ET';

% Cell edges
[XIe,ETe] = meshgrid(xievec,etevec);
XIe = XIe'; ETe = ETe';

% define our physical grid
XX = XI.*cos(ET);
YY = XI.*sin(ET);

% Transformation derivatives:
dx_dxi = @(xi,et)(      cos(et) );
dx_det = @(xi,et)( -xi.*sin(et) );
dy_dxi = @(xi,et)(      sin(et) );
dy_det = @(xi,et)(  xi.*cos(et) );

% Jacobian
Jt = XI;

% Save all the grid parameters (including edge values and tranformation functions)
grid_params.mxi = mxi; grid_params.met = met;
grid_params.dxi = dxi; grid_params.det = det;
grid_params.XX = XX;  grid_params.YY = YY;
grid_params.XI = XI;  grid_params.ET = ET;
grid_params.XIe = XIe; grid_params.ETe = ETe;
grid_params.dx_dxi = dx_dxi;
grid_params.dx_det = dx_det;
grid_params.dy_dxi = dy_dxi;
grid_params.dy_det = dy_det;

% Initial conditions
% qinit = .1 + sin(4*pi.*XX) + cos(2*pi.*YY);
xc = -0.5; yc = 0.0;
Rc = sqrt((XX-xc).^2 + (YY-yc).^2);
qinit = zeros(size(XX));
for i = 1:mxi
    for j = 1:met
        if (Rc(i,j)< 0.3)
            qinit(i,j) = cos(5/3*pi*Rc(i,j))^6;
        end
    end
end
q = qinit;

% Transformed variables
qt = q.*Jt;

% ---- Main time stepping loop ---- %
t = 0;  dT = Tf/nframes;
plot_stuff; pause(0.1);
for nf = 1:nframes

    qt = curvisolve_scalar(qt,t,dT,cfl_des,grid_params,false);
    t  = t+dT;
    q  = qt./Jt;

    % plot stuff (include the pause otherwise second figure doesn't always display)
    plot_stuff;  pause(0.1);

end
