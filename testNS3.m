clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set:
% (0) decide whether restar (restart=1) or scratch (restart=0)
% (1) the number of points (number of pressure cells)
% (2) the dimension of the domain
% (3) the Reynolds number and other geometrical characteristics
% (4) (for the moment) a fixed time step
%      and the number of time steps
% (5) the boundary values and rkpm yes/no
% (6) check convergence every iconv time steps
%     eventually save the field and
%     decide whether or not open "on the fly" figures      
%                                                                   
%               Ub(3) and Vb(3)                                      
%           ¡---------------------¡                                   
%   Ub(4)   ¡                     ¡ Ub(2) and Vb(2)                    
%   Vb(4)   ¡                     ¡                                     
%           -----------------------                                      
%               Ub(1) and Vb(1)                                           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (0)
restart=0
% (1)
ncx=850;
ncy=700;
% (2)
lx=17; 
ly=14;
% (3)
Re=30;
D=1;   %Diameter of the cylinder;
rc=D/2;
xc=3; % center of circle
yc=7;
% (4)
nt=20000; CFL=0.02;
% decide if rkpm or not
%irkpm=0 % no rkpm
irkpm=1 % yes rkpm
% (6) check conv
iconv=20;
snapshot_it=100;
nsteps=[1:snapshot_it:nt];
yes_snap=0;
[Ub,Vb]=set_velocity_bc(ncx,ncy,lx,ly);
% assemble the grid and the initial solution
if restart==0
   [xu,xv,xp,yu,yv,yp,u,v,p]=initialize(ncx,ncy,lx,ly,Ub,Vb);
   save xp.mat xp;
   save yp.mat yp;
   save xu.mat xu;
   save yu.mat yu;
   save xv.mat xv;
   save yv.mat yv;
else
   load uu2; load vv2; load pp2;
   load xp; load yp; load xu; load xv; load yu; load yv;
end
% initialize the matrices for Helmholz and Poisson Problems
nxu=length(xu); nyu=length(yu);
nxv=length(xv); nyv=length(yv);
% no uniform mesh with lx assumed not to be equal to ly!!!!
hx=lx/ncx; hy=ly/ncy;
Deltat=CFL*hx/max(Ub(4,:));
disp('Delta t ')
disp(Deltat)
beta=Deltat/hx;
alpha=Re/Deltat; 
%h=xp(2)-xp(1);
tic
[PUX,PUY,EUX,EUY]=assmatUFD(nxu,nyu,hx,hy);
[PVX,PVY,EVX,EVY]=assmatVFD(nxv,nyv,hx,hy);
[PX,PY,PMIX]=assmatFD(ncx,ncy,hx,hy);
toc
disp('finished eig problems ')
[xlag,ylag,theta,dS]=initialize_lag(xu,yu,rc,xc,yc,irkpm,0);
if irkpm > 0
tic 
[epsilonx,bcoefx,ix_x,jy_x,area_x,nsup_x,hx_x,hy_x]=precompute_eps(xu(2:length(xu)-1),yu(2:length(yu)-1),xlag,ylag,dS);
%keyboard
[epsilony,bcoefy,ix_y,jy_y,area_y,nsup_y,hx_y,hy_y]=precompute_eps(xv(2:length(xv)-1),yv(2:length(yv)-1),xlag,ylag,dS);
toc
else
[epsilonx,bcoefx,ix_x,jy_x,area_x,nsup_x,hx_x,hy_x]=precompute_eps_norkpm(xu(2:length(xu)-1),yu(2:length(yu)-1),xlag,ylag,dS);
[epsilony,bcoefy,ix_y,jy_y,area_y,nsup_y,hx_y,hy_y]=precompute_eps_norkpm(xv(2:length(xv)-1),yv(2:length(yv)-1),xlag,ylag,dS);
end
disp('finished precomputing immersed boundary data')
%keyboard
%%%% start the time loop
disp('starting time advancement')
isnapshot=0;
Qin=sum(u(nxu,2:nyu-1)); uslice=zeros(2,nyu); vslice=zeros(2,nyv);
for it=1:nt
disp('step number'); disp(it)
tic
    actime=(it-1)*Deltat;
    if mod(it,iconv) == 0
       uold=u; vold=v;
    end
%%%%%%%% X_MOM EQUATION %%%%%%%%%%%%%%
% compute the laplacian of u
%   lapla=laplacian(u,hx,hy,Ub,1);
% compute NL terms of X-mom equation
   NL=nonlinear(u,v,Ub,Vb,lx,ly,ncx,ncy,1);
% set the rhs for the Helm problem
   rhs=Re*NL-alpha*u(2:nxu-1,2:nyu-1);
   [Ubs]=setBCout(u,Deltat,hx,Qin,1);
% solve the Helm problem
   ustar=solveUFD(PUX,PUY,EUX,EUY,rhs,Ub,Ubs,hx,hy,alpha);
% compute and spread forces
   forcex=0;
   rhs=immersed_boundary(-ustar/Deltat,xu(2:length(xu)-1),yu(2:length(yu)-1),xlag,ylag,dS,epsilonx,bcoefx,...
   ix_x,jy_x,area_x,nsup_x,hx_x,hy_x); 
   forcex=forcex+hx*hy*sum(sum(rhs));
   U1=zeros(4,max(nxu,nyu)); U2=zeros(nyu,1);
   rhs=solveUFD(PUX,PUY,EUX,EUY,-Re*rhs,U1,U2,hx,hy,alpha);
   ustar=ustar+rhs;
   Cd=2*forcex/D;
%   correct the eventual non uniformity
%%%%%%%% Y_MOM EQUATION %%%%%%%%%%%%%%
% compute the laplacian of v
%   lapla=laplacian(v,hx,hy,Vb,2);
% compute NL terms of Y-mom equation
   NL=nonlinear(u,v,Ub,Vb,lx,ly,ncx,ncy,2);
   rhs=Re*NL-alpha*v(2:nxv-1,2:nyv-1);
   vstar=solveVFD(PVX,PVY,EVX,EVY,rhs,Vb,hx,hy,alpha);
% compute the force for the y-mom equation
   forcey=0;
   rhs=immersed_boundary(-vstar/Deltat,xv(2:length(xv)-1),yv(2:length(yv)-1),xlag,ylag,dS,epsilony,bcoefy,...
   ix_y,jy_y,area_y,nsup_y,hx_y,hy_y); 
      forcey=forcey+hx*hy*sum(sum(rhs));
   V1=zeros(4,max(nxv,nyv));
   rhs=solveVFD(PVX,PVY,EVX,EVY,-Re*rhs,V1,hx,hy,alpha);
   vstar=vstar+rhs;
   Cf=2*forcey/D;
% prepare the rhs for the pressure equation
   [u,v]=setbc(ustar,vstar,Ub,Vb,Ubs);
   div=divergence(u,v,hx,hy,ncx,ncy);
% solve the Poisson problem for pressure
     phi=solveFD(PX,PY,PMIX,div/Deltat);
%%%MC p=fastsolver_P(PP,PP1,PL,div);
% projection step
  [u,v,p]=finalvelocity(phi,u,v,Deltat,hx,hy);
%  keyboard
  if mod(it,iconv) == 0
     isnapshot=isnapshot+1;
     disp('/****** at time ***/');
     actime
     difu=abs(uold-u); difv=abs(vold-v);
     difumax=max(max(difu));
     difvmax=max(max(difv));
     divmax=divergence_check(u,v,hx,hy,ncx,ncy);
     vinlet=sum(u(1,2:nyu-1));
     voutlet=sum(u(nxu,2:nyu-1));
     disp('max div');
     disp(divmax)
     disp('max dif in u');
     disp(difumax)
     disp('max dif in v');
     disp(difvmax)
     disp('CD');
     disp(Cd)
     disp('CL');
     disp(Cf)
     disp('sum vel inlet');
     disp(vinlet)
     disp('sum vel outlet');
     disp(voutlet)
     disp('saving the field')
     
     save uu.mat u;
     save vv.mat v;
     save pp.mat p;
     forcexx(isnapshot)=forcex;
     forceyy(isnapshot)=forcey;
     Cdd(isnapshot)=Cd;
     Cll(isnapshot)=Cf;
     time(isnapshot)=actime;
     save time.mat time;
     save forcexx.mat forcexx;
     save forceyy.mat forceyy;
     save Cdd.mat Cdd;
     save Cll.mat Cll;
  end
toc

disp('time for 1 time step')
 % snapshots  
%     if find(nsteps-it==0)
%           fprintf('\t\tsnapshot:%d\n',it);
%    	figure(1);clf;hold on;axis equal;pcolor(xu,yu,u');plot(xlag,ylag,'o-r');
%       shading interp;caxis([-1 1]);colormap(jet2);colorbar;axis([0 lx 0 ly]);drawnow;
%   	figure(2); clf;hold on;axis equal;pcolor(xv,yv,v');plot(xlag,ylag,'o-r');
%    	shading interp;caxis([-1 1]);colormap(jet2);colorbar;axis([0 lx 0 ly]);drawnow;
%         % figure(3); clf;hold on;axis equal;pcolor(X,Y,omega(2:ncx,2:ncy)');plot(xlag,ylag,'o-r');
%    	%shading interp;caxis([-1 1]);colormap(jet2);colorbar;axis([0 lx 0 ly]);drawnow;
%   end
end
[omega,psi,xo,yo]=vort(u,v,hx,hy);
[X,Y]=meshgrid(xo(2:ncx),yo(2:ncy));
%[omega,psi,xo,yo]=vort(u,v,hx,hy);
fp=figure;
ap=axes;
maxval=max(max(psi));
minval=min(min(psi));
vals=linspace(minval,maxval,50);
contour(X,Y,psi',vals);
set(ap,'DataAspectRatio',[1 1 1],'XLim',[0 lx],'YLim',[0 ly]);
title(sprintf(' Stream function  (%g...%g); Re_b=%g',...
                minval,maxval,Re));
colorbar;

fo=figure;
ao=axes;
maxval=max(max(omega));
minval=min(min(omega));
vals=linspace(minval,maxval,100);
contour(X,Y,omega(2:ncx,2:ncy)',vals);
set(ao,'DataAspectRatio',[1 1 1],'XLim',[0 lx],'YLim',[0 ly]);
title(sprintf(' Vorticity (%g...%g); Re_b=%g',...
                minval,maxval,Re));
colorbar;

fx=figure;
ax=axes;
maxval=max(max(u));
minval=min(min(u));
vals=linspace(minval,maxval,50);
contour(xu,yu,u',vals);
set(ax,'DataAspectRatio',[1 1 1],'XLim',[0 lx],'YLim',[0 ly]);
title(sprintf(' u-com velocity (%g...%g); Re_b=%g',...
                minval,maxval,Re));
colorbar;


