clear all
iflap=0;
ierror=0;
%%%% INPUT SECTION %%%%%%%
ncx=500; ncy=420;
lx=16.5; ly=13.86;
xc=3.003; yc=6.930; rc=0.5;
D=2*rc; Re=100;
%%%% LOAD SECTION %%%%%%%
load uu; load vv; load pp; load phi;
load ulag; load vlag; load plag; load velnorm;
%%%% GRID SECTION %%%%%%
hx=lx/ncx; hy=ly/ncy;
nxu=ncx+1; nyu=ncy+2;
nxv=ncx+2; nyv=ncy+1;
for i=1:nxu
    xu(i)=(i-1)*hx;
end
for j=1:nyu
    yu(j)=-hy/2+(j-1)*hy;
end
for i=1:nxv
    xv(i)=-hx/2+(i-1)*hx;
end
for j=1:nyv
    yv(j)=(j-1)*hy;
end
xp=xv(2:ncx+1); yp=yu(2:ncy+1);
%[xlag,ylag,theta,dS]=initialize_lag(xu,yu,rc,xc,yc,1,0);
load ango.mat;
laango=length(ango);
%[xlag,ylag,theta,dS,nss,nse]=initialize_lag_tail(xu,yu,rc,xc,yc,0,0,ango(laango));
[xlag,ylag,theta,dS,nss,nse]=initialize_lag_tail(xu,yu,rc,xc,yc,0,0,0);
%%%%% ERROR SECTION  %%%%
if ierror==1
veltan=-ulag.*sin(theta)'+vlag.*cos(theta)';
[vnmax,inm]=max(abs(velnorm)); tinm=radtodeg(theta(inm));
[vtmax,itm]=max(abs(veltan));  titm=radtodeg(theta(itm));
[divmax,idm]=max(abs(plag));  tidm=radtodeg(theta(idm));
fprintf(1,'Max norm vel=%8.3e ; Location max norm (degrees) =%8.3f \n',vnmax,tinm);
fprintf(1,'Max tang vel=%8.3e ; Location max tang (degrees) =%8.3f \n',vtmax,titm);
fprintf(1,'Max divergence on Gamma=%8.3e ; Location max divergence (degrees) =%8.3f \n',divmax,tidm);
%%%%%%%%%%%%%%%
f1=figure; %a1=axes;
vlim2=1.25*max(velnorm); vlim1=1.25*min(velnorm);
plot(radtodeg(theta),velnorm);
axis([0 365 vlim1 vlim2]);
set(gca,'XTick',[0 45 90 135 180 225 270 315 360])
%%%%%%%%%%%%%%
f2=figure; %a1=axes;
vlim2=1.25*max(veltan); vlim1=1.25*min(veltan);
plot(radtodeg(theta),veltan);
axis([0 365 vlim1 vlim2]);
set(gca,'XTick',[0 45 90 135 180 225 270 315 360])
%%%%%%%%%%%%%%
f3=figure; %a1=axes;
vlim2=1.25*max(plag); vlim1=1.25*min(plag);
plot(radtodeg(theta),plag);
axis([0 365 vlim1 vlim2]);
set(gca,'XTick',[0 45 90 135 180 225 270 315 360])
end
%%%%%%%%% PLOT SECTION
ic=ceil(xc/hx)+1; jc=ceil(yc/hy)+1;
ir1=ceil((xc-rc)/hx); ir2=ceil((xc+rc)/hx);
idi=ir2-ir1;
istart=ic-floor(0.75*idi); iend=ic+floor(7.75*idi);
jstart=jc-floor(1.5*idi); jend=jc+floor(1.5*idi);
%
f4=figure; a1=axes;
[X,Y]=meshgrid(xu(istart:iend),yu(jstart:jend));
minval=min(min(u)); maxval=max(max(u));
conts=linspace(minval,maxval,100);
contour(X,Y,u(istart:iend,jstart:jend)',conts);
set(a1,'DataAspectRatio',[1 1 1],'XLim',[xu(istart) xu(iend)],'YLim',[yu(jstart) yu(jend)]);
axis equal
title(sprintf(' Stream wise velocity  (%g...%g); Re_b=%g',...
                minval,maxval,Re));
colorbar;
hold on; plot(xlag,ylag,'k.');
%
f5=figure; a2=axes;
[X,Y]=meshgrid(xv(istart:iend),yv(jstart:jend));
minval=min(min(v)); maxval=max(max(v));
conts=linspace(minval,maxval,100);
contour(X,Y,v(istart:iend,jstart:jend)',conts);
set(a2,'DataAspectRatio',[1 1 1],'XLim',[xv(istart) xv(iend)],'YLim',[yv(jstart) yv(jend)]);
axis equal
title(sprintf(' Normal velocity  (%g...%g); Re_b=%g',...
                minval,maxval,Re));
colorbar;
hold on; plot(xlag,ylag,'k.');
%
f6=figure; a3=axes;
[X,Y]=meshgrid(xp(istart:iend),yp(jstart:jend));
minval=min(min(p)); maxval=max(max(p));
conts=linspace(minval,maxval,100);
contour(X,Y,p(istart:iend,jstart:jend)',conts);
set(a3,'DataAspectRatio',[1 1 1],'XLim',[xp(istart) xp(iend)],'YLim',[yp(jstart) yp(jend)]);
axis equal
title(sprintf(' Pressure (%g...%g); Re_b=%g',...
                minval,maxval,Re));
colorbar;
hold on; plot(xlag,ylag,'k.');
%
if ierror==1
  f7=figure; a4=axes;
  [X,Y]=meshgrid(xp(istart:iend),yp(jstart:jend));
  minval=min(min(phi)); maxval=max(max(phi));
  conts=linspace(minval,maxval,100);
  contour(X,Y,phi(istart:iend,jstart:jend)',conts);
  set(a4,'DataAspectRatio',[1 1 1],'XLim',[xp(istart) xp(iend)],'YLim',[yp(jstart) yp(jend)]);
  axis equal
  title(sprintf(' Pressure Correction (%g...%g); Re_b=%g',...
                minval,maxval,Re));
  colorbar;
  hold on; plot(xlag,ylag,'k.');
%
[n,m]=size(phi);
  mask=ones(n,m);
  for i=1:n
     d1=(xp(i)-xc)^2;
     for j=1:m
         d2=(yp(j)-yc)^2;
         if sqrt(d1+d2) > rc
            mask(i,j)=0;
         end
     end
  end
%
  f8=figure; a5=axes;
  [X,Y]=meshgrid(xp(istart:iend),yp(jstart:jend));
  minval=min(min(mask.*phi)); maxval=max(max(mask.*phi));
  conts=linspace(minval,maxval,100);
  contour(X,Y,(mask(istart:iend,jstart:jend).*phi(istart:iend,jstart:jend))',conts);
  set(a5,'DataAspectRatio',[1 1 1],'XLim',[xp(istart) xp(iend)],'YLim',[yp(jstart) yp(jend)]);
  axis equal
  title(sprintf(' Pressure Correction inside body (%g...%g); Re_b=%g',...
                minval,maxval,Re));
  colorbar;
end
%%%%%%%%%%%
f9=figure; a6=axes;
[X,Y]=meshgrid(xp(istart:iend),yp(jstart:jend));
[omega,psi,xo,yo]=vort(u,v,hx,hy);
minval=min(min(omega)); maxval=max(max(omega));
lung=0.8*(maxval-minval); inc=lung/400;
conts=0.8*minval:4*inc:0.8*maxval;
%conts=linspace(minval,maxval,100);
contour(X,Y,(omega(istart:iend,jstart:jend))',conts);
set(a6,'DataAspectRatio',[1 1 1],'XLim',[xp(istart) xp(iend)],'YLim',[yp(jstart) yp(jend)]);
axis equal
title(sprintf(' vorticity  (%g...%g); Re_b=%g',...
                minval,maxval,Re));
colorbar;
hold on; plot(xlag,ylag,'k.');
%%%%%%%%%%%%%%
if iflap==1
  load ango; load ango_dot; load torques; load time;
  f10=figure; plot(time,ango); title('angle vs time');
  f11=figure; plot(time,ango_dot); title('angular velocity vs time');
  f12=figure; plot(time,torques); title('fluid torque vs time');
end
