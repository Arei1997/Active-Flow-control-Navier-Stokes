clear all
load pp;
load xp;
load yp;
ncx=length(xp); 
ncy=length(yp);
%%%%%%% INPUT SECTION
lx=17; % streamwise domain extension
ly=14; % spanwise domain extension 
xc=3;  % x-center of circle
yc=7;  % y-center of circle
Re=30; % Reynolds number
D=1;  %Diameter of the cylinder;
%%%%%%% END INPUT SECTION
rc=D/2; 
dxmin=min(diff(xp));
dymin=min(diff(yp));
diagm=sqrt(dxmin^2+dymin^2);
nls=floor(0.95*2*pi*rc/diagm);
[theta,xlag,ylag,dS]=setupdisklag(nls,xc,yc,rc,0);
xx=xp(2:length(xp)-1); yy=yp(2:length(yp)-1);
[epsp,bcoefp,ixp,jyp,areap,nsupp,hxp,hyp]=precompute_eps(xx,yy,xlag,ylag,dS);
%[epsp,bcoefp,ixp,jyp,areap,nsupp,hxp,hyp]=precompute_eps_norkpm(xx,yy,xlag,ylag,dS);
plag=interpol(p,nls,ncx-2,ncy-2,xlag,ylag,dS,xx,yy,bcoefp,ixp,jyp,areap,nsupp,hxp,hyp);
nstart=ceil(nls/2);
plot(theta(nstart:nls),abs(plag(nstart:nls)),'r',theta(nstart-1:-1:1),abs(plag(nstart-1:-1:1)),'b');
plim=1.25*max(abs(plag));
axis([0 2*pi 0 plim]);
legend('bottom half', 'top half')
title('\theta=0 is the trailing edge position')
%%%%%%%%%%%%%%%%%%%%%%
load uu; load vv;
hx=xp(2)-xp(1); hy=yp(2)-yp(1);
[omega,psi,xo,yo]=vort(u,v,hx,hy);
[nx1,ny1]=size(psi);
for i=1:nx1
    psi(i,:)=psi(i,:)+yo(2:length(yo)-1);
end
psilag=interpol(psi,nls,ncx-2,ncy-2,xlag,ylag,dS,xx,yy,bcoefp,ixp,jyp,areap,nsupp,hxp,hyp);
figure
plot(theta(nstart:nls),abs(psilag(nstart:nls)),'r',theta(nstart-1:-1:1),abs(psilag(nstart-1:-1:1)),'b');
plim=1.25*max(abs(psilag));
axis([0 2*pi 0 plim]);
legend('bottom half', 'top half')
title('\theta=0 is the trailing edge position')
