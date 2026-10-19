function [xlag,ylag,epx,epy,bx,by,ix1,ix2,jy1,jy2,a1,a2,nsx,nsy,hx1,hx2,hy1,hy2]=...
update_lag_onlytail(xu,yu,xv,yv,dS,nls_s,nls_e,xlag,ylag,xc,yc,crot,alpha)
% lagrangian parameters
dxmin=min(diff(xu));
dymin=min(diff(yu));
diagm=sqrt(dxmin^2+dymin^2);
%
lbeam=0;
for i=nls_s:nls_e-1
   lbeam=lbeam+dS(i);
end
nls=nls_e-nls_s+1;
nodes=linspace(0,lbeam,nls);
for i=nls_s:nls_e
    xlag(i)=xc-crot*cos(alpha)+nodes(i)*cos(alpha);
    ylag(i)=yc-crot*sin(alpha)+nodes(i)*sin(alpha);
end
[epx,bx,ix1,jy1,a1,nsx,hx1,hy1]=precompute_eps_norkpm(xu(2:length(xu)-1),yu(2:length(yu)-1),xlag,ylag,dS);
[epy,by,ix2,jy2,a2,nsy,hx2,hy2]=precompute_eps_norkpm(xv(2:length(xv)-1),yv(2:length(yv)-1),xlag,ylag,dS);
%[epx,bx,ix1,jy1,a1,nsx,hx1,hy1]=precompute_eps(xu(2:length(xu)-1),yu(2:length(yu)-1),xlag,ylag,dS);
%[epy,by,ix2,jy2,a2,nsy,hx2,hy2]=precompute_eps(xv(2:length(xv)-1),yv(2:length(yv)-1),xlag,ylag,dS);
