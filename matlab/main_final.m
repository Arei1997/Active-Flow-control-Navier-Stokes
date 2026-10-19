clear all
global nls ni nj lx ly xc yc rc
global xlag ylag dS flag x y b stop

%  -  Solve heat equation around a fixed/moving cylinder
%  -  finite difference in space and euler implicit in time
%  -  using immersed boundaries imposed by RKPM
%  -  compute the correct deltav for the spreading step
%  -  17 november 2009

!rm -f check_interp.dat
ni=71; nj=71; ntot=(ni-2)*(nj-2);
tmax1=1; tmax2=0; tmax3=0; tmax4=0;  tmax=tmax1;     %only left
%tmax1=3; tmax2=60; tmax3=12; tmax4=6;  tmax=tmax1;  %everywhere
v=linspace(-tmax/10,tmax/2,20);
lx=1; ly=1; facx=3; facy=2;
x=nu_grid(ni,lx,facx,0);
y=nu_grid(nj,ly,facy,0);
[A]=genmat(x,y);
[X,Y]=meshgrid(x(2:ni-1),y(2:nj-1));
deltat=5e-6; kappa=5.2; nmax=200;
%%%%%%%%%%%%%% SETUP SOLID BODY
rc=0.25; xc=0.32; yc=0.5; nls=115;
irkpm=3;
%irkpm=conserve 0th integral (1) linear moments (2) 2nd-order moment (3), Roma (0)
[theta]=setupdisk_fourier();
A=A+speye(size(A))/kappa/deltat;
T=zeros(length(A(:,1)),1); q=zeros(length(A(:,1)),1);

mindx=min(min(diff(x)),min(diff(y)));
minds=2*pi*rc/nls;
fprintf('min dx:%d, min ds:%d\n',mindx,minds);
figure(798);clf;pcolor(X,Y,NaN(size(X)));hold on; plot(xlag,ylag,'o-r'); shading faceted;axis equal

%%%%%%%%%%%%% BOUNDARY CONDITIONS
    ii=2;    % WEST , x=0
    dxm=x(ii)-x(ii-1);dxp=x(ii+1)-x(ii);
    factorSX=2/(dxp+dxm)/dxm;
    q(1:nj-2)=tmax1*sin(pi*y(2:nj-1))*factorSX;
    ii=ni-1; % EAST, x=L
    dxm=x(ii)-x(ii-1);dxp=x(ii+1)-x(ii);
    factorDX=2/(dxp+dxm)/dxp;
    q(end-(nj-2)+1:end)=tmax2*sin(pi*y(2:nj-1))*factorDX;
    jj=2;    % BOTTOM, y=0
    dxm=y(jj)-y(jj-1); dxp=y(jj+1)-y(jj);
    factorBOT=2/(dxp+dxm)/dxm;
    q(1:(nj-2):end)=tmax3*sin(pi*x(2:ni-1))*factorBOT;
    jj=nj-1; % TOP, y=1
    dxm=y(jj)-y(jj-1);dxp=y(jj+1)-y(jj);
    factorTOP=2/(dxp+dxm)/dxp;
    q((nj-2):(nj-2):end)=tmax4*sin(pi*x(2:ni-1))*factorTOP;
%%%%%%%%%%%%%%% LOOP IN TIME IMPLICIT EULER
epsilon(1:nls)=0;
for n=1:nmax
    disp(n)
    rhs=T/kappa/deltat+q;
    T=A\rhs;
    % steps in physical formalism
    fun=zeros(nj,ni);
    fun(2:nj-1,2:ni-1)=reshape(T,nj-2,ni-2); fun=fun/deltat;

    % uncomment for moving cylinder
%     xlag=xlag ;%+ 0.01*cos(pi*n/15); 
%     ylag=ylag+ 0.01*sin(pi*n/15);

%     b=assemble_rkpm();
     b=zeros(6,nls);b(1,:)=1;

    % interpol eulerian field on the lagrangian points
    [val]=interpol(fun);

% FIND EPSILON
%-------------------------------
[flag]=ones(nls,1);
%[flag]=interpol(u);
bic_b=interpol(fun);
%------------------------------- matrix formulation
[M,testdiag]=buildmat;%figure(244);clf;imagesc(M);colorbar % en dur
eta=M\ones(nls,1);
epsilon=eta;
%-------------------------------

tolit=min([min(diff(x)),min(diff(y))])^2;
%guess=ones(1,nls);
guess=zeros(1,nls);
maxit=nls;
%bic_b=ones(nls,1);
restart_gmres=nls;
precond_bicg=1;
gmres_params=[tolit,maxit];

%fprintf('iterative method...\n')
%[eta,flag2,relres2,iter2,resvec2]=gmres(@matvec,bic_b,restart_gmres,tolit,maxit,@precondfun,[],guess(:));

%eta=M\bic_b;
%diag_mat=precondfun(bic_b);
%for ix=1:nls
%	if abs(bic_b(ix))>1e-6/deltat
%	epsilon(ix)=1/bic_b(ix)*eta(ix);
%        else
%%		disp 'ok'
%%	epsilon(ix)=diag_mat(ix);
%	epsilon(ix)=guess(ix);
%	end
%end	
%epsilon=epsilon(:);

figure(245),clf
plot(epsilon);drawnow

    % spread the lagrangian forcing on the eulerian field
    [force]=spread(-val,epsilon);

%    % uncomment to check the quality of interpolation
%     [valtest]=interpol(-force);
%     figure(2);clf;plot(val,'ob');hold on; plot(valtest,'r')
%     fprintf('error on lagragian forcing: %1.5e\n',norm(val-valtest))
%     pause

    % steps in algebraic formalism
    ff=reshape(force(2:nj-1,2:ni-1),ntot,1);
    rhs=rhs+ff/kappa;
    T=A\rhs;
    % steps in physical formalism
    t=reshape(T,nj-2,ni-2);
    figure(10);contourf(X,Y,t,v); shading faceted; axis equal; hold on; plot(xlag,ylag,'o-r');
    caxis([-1 1]);colormap(jet2);
%     pcolor(X,Y,t); axis equal; hold on; plot(xlag,ylag,'w*');pause(0.01);hold off;shading faceted
    drawnow;pause(0.01);hold off

    % uncomment to save snaphots
    % nom=sprintf('ANIM/%d.png',n);
    % print('-dpng',nom);

end

