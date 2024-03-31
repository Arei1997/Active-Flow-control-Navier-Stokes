function CD = solver(angle)

ncx = 500; ncy = 420;     % number of nodes for the pressure in x and y

lx = 16.5; ly = 13.86;    % x-wise and y-wise extension of the domain 

Re = 150;                 % Reynolds number of the immersed circle (U upstream * Diameter/kinematic visc.)
D = 1; rc = D/2;          % Diameter of the cylinder, circle radius
xc = 3.003; yc = 6.930;   % center of circle (coordinates)

nt = 50; CFL = 0.2;     % number of time steps and CFL number

[Ub, Vb] = set_velocity_bc(ncx, ncy, lx, ly);               % set velocity on the boundaries

% set the initial flap position and constants
length_beam = D; 
center_rot = length_beam; 
k_spring = 0.5; 
rho_beam = 2; Delta_rho = rho_beam - 1;
Inertia_mom = Delta_rho*length_beam^3/3;
Inertia_mom = Inertia_mom + Delta_rho*length_beam*center_rot*(center_rot - length_beam);

% assemble the grid and the initial solution
angolo = angle; angolo_dot = 0; TORQUE = 0;
[xu, xv, xp, yu, yv, yp, u, v, p] = initialize(ncx, ncy, lx, ly, Ub, Vb); % initialize the field and grid if from scratch

% initialize the matrices for Helmholz and Poisson Problems
nxu = length(xu); nyu = length(yu);            % no. of nodes for u (x-comp. of velocity)
nxv = length(xv); nyv = length(yv);            % no. of nodes for v (y-comp. of velocity)
nxp = length(xp); nyp = length(yp);            % no. of nodes for p 

% no uniform mesh with lx assumed not to be equal to ly!!!!
hx = lx/ncx; hy = ly/ncy;                      % mesh spacing
Deltat = CFL*hx/max(Ub(4,:));                  % Allowed time step using entrance velocity as max
beta = Deltat/hx; alpha = Re/Deltat;           % some parameters 

% D'Yakonov scheme (approximate factorization of viscous terms)
[alxu, bexu, cxu, alyu, beyu, cyu] = predyakonovU(nxu, nyu, hx, hy, Deltat, Re);
[alxv, bexv, cxv, alyv, beyv, cyv] = predyakonovV(nxv, nyv, hx, hy, Deltat, Re);
[PX, PY, PMIX] = assmatFD(ncx, ncy, hx, hy);    % diagonalize pressure laplacian operator

[xlag, ylag, theta, dS, nlsts, nlste] = initialize_lag_tail(xu, yu, rc, xc, yc, 0, 0, angolo); 

% compute the delta mollifier and the width of the corona around the cirlce (\epsilon in Jou. Comp. Phys. Paper)
[epsilonx, bcoefx, ix_x, jy_x, area_x, nsup_x, hx_x, hy_x] = precompute_eps_norkpm(xu(2:length(xu)-1), yu(2:length(yu)-1), xlag, ylag, dS);
[epsilony, bcoefy, ix_y, jy_y, area_y, nsup_y, hx_y, hy_y] = precompute_eps_norkpm(xv(2:length(xv)-1), yv(2:length(yv)-1), xlag, ylag, dS);

% start the time loop
isnapshot = 0;
Qin = sum(u(nxu, 2:nyu-1));                      % flow rate at inlet
AB1 = 1; AB2 = 0; CGP = 0;                       % first step Euler for nonlinear terms, from the second one Adams-Bashfort
laplax = zeros(nxu-2, nyu-2); NLx = zeros(nxu-2, nyu-2);
laplay = zeros(nxv-2, nyv-2); NLy = zeros(nxv-2, nyv-2);
icount = 0;

for it = 1:nt
    actime = (it-1)*Deltat;
    if mod(it, 20) == 0
        uold = u; vold = v;
    end

    %%% X_MOM EQUATION %%%
    % FROM A to B prediction of u from x momentum equation without immersed object (fully explicit)
    
    % --------- A --------------
    lapla = laplacian(u, hx, hy);                % compute the laplacian of u
    NL = nonlinear(u, v, Ub, Vb, lx, ly, ncx, ncy, 1);                       % NL terms of X-mom equation 
    GP = pressgrad(p, hx, hy, 1);                                            % dp/dx
    
    % compute the force for the x-mom equation
    rhs = sumupAB(u, lapla, NL, laplax, NLx, CGP*GP, Deltat, Re, AB1, AB2);  % (2.10a)
    velx = -angolo_dot*(ylag(nlsts:nlste)-yc);
    [rhs, val_x] = immersed_boundary(rhs, xu(2:length(xu)-1), yu(2:length(yu)-1), xlag, ylag,...
            dS, epsilonx, bcoefx, ix_x, jy_x, area_x, nsup_x, hx_x, hy_x, nlsts, nlste, velx/Deltat);        
    FX = -val_x(nlsts:nlste).*epsilonx(nlsts:nlste).*dS(nlsts:nlste)';
    
    % --------- B -------------- 
    % rhs contains the spread force field required to impose bc's on the circle for u
    forcex = hx*hy*sum(sum(rhs)); Cd = 2*forcex/D;     % compute integral force in x direction on the circle and the drag coeficient
    % assemble rhs for the implicit step
    % from C to D redo x-momentum with the spreaded body force on the rhs 
    % Approximate fact. for viscous terms and Adams Bashfort for the non-linear ones
    
    %---------C ------------
    rhs = rhs - AB1*NL + AB2*NLx - CGP*GP;
    [Ubs] = setBCout(u, v, Deltat, hx, Qin, 1);
    ustar = DyakonovU(rhs, u, Ub, Ubs, alxu, bexu, cxu, alyu, beyu, cyu, hx, hy, Re, Deltat);
    
    %---------D ------------
    laplax = lapla; NLx = NL;
    % DO STAGE A-B and C-D for the y momentum equation
    
    %%% Y_MOM EQUATION %%%
    lapla = laplacian(v, hx, hy);                 % compute the laplacian of v
    NL = nonlinear(u, v, Ub, Vb, lx, ly, ncx, ncy, 2); % NL terms of Y-mom equation 
    GP=pressgrad(p,hx,hy,2);                            % dp/dy
    
    % compute the force for the y-mom equation
    rhs = sumupAB(v, lapla, NL, laplay, NLy, CGP*GP, Deltat, Re, AB1, AB2);
    vely = angolo_dot*(xlag(nlsts:nlste) - xc);
    [rhs, val_y] = immersed_boundary(rhs, xv(2:length(xv)-1), yv(2:length(yv)-1), xlag, ylag,...
           dS, epsilony, bcoefy, ix_y, jy_y, area_y, nsup_y, hx_y, hy_y, nlsts, nlste, vely/Deltat); 
    FY = -val_y(nlsts:nlste).*epsilony(nlsts:nlste).*dS(nlsts:nlste)';
    forcey = hx*hy*sum(sum(rhs)); Cf = 2*forcey/D;
    
    % assemble rhs for the implicit step
    rhs = rhs - AB1*NL + AB2*NLy - CGP*GP;
    [Vbs] = setBCout(u, v, Deltat, hx, Qin, 2);
    vstar = DyakonovV(rhs, v, Vb, Vbs, alxv, bexv, cxv, alyv, beyv, cyv, hx, hy, Re, Deltat);
    laplay = lapla; NLy = NL;
    
    %%% Get ready for the projection step %%%
    [u, v] = setbc(ustar, vstar, Ub, Vb, Ubs, Vbs);      % set the boundary values on the boundary frame
    div = divergence(u, v, hx, hy, ncx, ncy);            % prepare the rhs for the pressure equation
    phi = solveFD(PX, PY, PMIX, div/Deltat);             % solve the Poisson problem for pressure
    [u, v, p] = finalvelocityVK(phi, u, v, p, Deltat, Re, hx, hy);
    ulag = interpolante(u(2:nxu-1,2:nyu-1), xu(2:length(xu)-1), yu(2:length(yu)-1), xlag, ylag,...
        dS, bcoefx, ix_x, jy_x, area_x, nsup_x, hx_x, hy_x);
    vlag = interpolante(v(2:nxv-1,2:nyv-1), xv(2:length(xv)-1), yv(2:length(yv)-1), xlag, ylag,...
        dS, bcoefy, ix_y, jy_y, area_y, nsup_y, hx_y, hy_y);
    v1(1) = max(abs(ulag)); v1(2) = max(abs(vlag)); 
    if mod(it, 20) == 0
        icount = icount + 1;        
        isnapshot = isnapshot + 1;
        dummy = screenout(uold, u, vold, v, actime, hx, hy, ncx, ncy, Cd, Cf);
        forcexx(isnapshot) = forcex; forceyy(isnapshot) = forcey;
        Cdd(isnapshot) = Cd; Cll(isnapshot) = Cf;
        time(isnapshot) = actime;
        dummy = savefield(u, v, p, phi, time, forcexx, forceyy, Cdd, Cll);
    end    
    % update the lagrangian nodes and related variables
    TORQUE = sum((xlag(nlsts:nlste)-xc).*FY'-(ylag(nlsts:nlste)-yc).*FX');
    AB1 = 3/2; AB2 = 1/2; CGP = 1; % Adam Bashfort coef.
end % END of time looping

CD =  Cdd(end);

