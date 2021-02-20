% Sample 2D advection solver on curvilinear coordinates

clear;

% define our curvilinear coordinates (or computational graph) values
mxi = 100; met = 101;                       % spacing

xi_min = 0; xi_max = 1;
et_min = 0; et_max = 1;

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

% perturbation values
ex = 0.01; ey = 0.02; ax = 2; ay = 4; 

% define our physical grid
XX = XI + ex*sin(2*pi.*ET.*ax);
YY = ET + ey*sin(2*pi.*XI.*ay);

% Transformation derivatives:
dx_dxi = @(xi,et)( 1.0*ones(size(xi)) );
dx_det = @(xi,et)( (2.0*pi*ax)*ex*cos(2.0*pi*ax*et) );
dy_dxi = @(xi,et)( (2.0*pi*ay)*ey*cos(2.0*pi*ay*xi) );
dy_det = @(xi,et)( 1.0*ones(size(xi)) );

% Save all the grid parameters (including edge values and tranformation functions)
grid.mxi = mxi; grid.met = met;
grid.dxi = dxi; grid.det = det;
grid.XX = XX;  grid.YY = YY;
grid.XI = XI;  grid.ET = ET;
grid.XIe = XIe; grid.ETe = ETe;
grid.dx_dxi = dx_dxi;
grid.dx_det = dx_det;
grid.dy_dxi = dy_dxi;
grid.dy_det = dy_det;

% problem specific parameters, in the physical domain
u = 2.0;      v = 1.0;

% problem specific parameters
params.u = u;
params.v = v;

% Initial conditions
% qinit = .1 + sin(4*pi.*XX) + cos(2*pi.*YY);
xc = 0.5; yc = 0.5;
Rc = sqrt((XX-xc).^2 + (YY-yc).^2);
qinit = zeros(size(XI));
for i = 1:mxi
    for j = 1:met
        if (Rc(i,j)< 0.3)
            qinit(i,j) = cos(5/3*pi*Rc(i,j))^6;
        end
    end
end
q = qinit;

% define our local Jacobian value Jt and qt (qt representing q tilda)
% TODO - need to compute what this tranformation actually is here
Jt = 1.0;
qt = q./Jt;

% plot the initial conditions in the physical and computational domain
figure(1);
clf
pcolor(XX,YY,q);
%contourf(XX',YY',q');
t1 = title(['Advection equation curvilinear example t = ' num2str(0)]);
set(t1,'FontSize',12);
xlabel('x', 'FontSize', 16); ylabel('y', 'FontSize', 16);

figure(2);
clf
pcolor(XI',ET',qt');
%contourf(XI',ET',qt');
t1 = title(['Advection equation curvilinear example - computational domain t = ' num2str(0)]);
set(t1,'FontSize',12);
xlabel('\xi','FontSize', 16); ylabel('\eta', 'FontSize', 16);
% pause;

% ---- Main time stepping loop ---- %
t = 0;  dT = Tf/nframes;
for nf = 1:nframes
    
    qt = curvisolve_scalar(qt,t,dT,cfl_des,grid,params);
    t  = t+dT;
    q  = qt.*Jt;
    
    figure(1);
    clf
    pcolor(XX',YY',q');
    %contourf(XX',YY',q');    
    t1 = title(['Advection equation curvilinear example t = ' num2str(nf*(Tf/nframes))]);
    set(t1,'FontSize',12);
    xlabel('x', 'FontSize', 16); ylabel('y', 'FontSize', 16);

    figure(2);
    clf
    pcolor(XI',ET',qt');
    %contourf(XI',ET',qt');    
    t1 = title(['Advection equation curvilinear example - computational domain t = ' num2str(nf*(Tf/nframes))]);
    set(t1,'FontSize',12);
    xlabel('\xi','FontSize', 16); ylabel('\eta', 'FontSize', 16);
    %    pause;
    
end
