function [xlag,ylag,epx,epy,bx,by,ix1,ix2,jy1,jy2,a1,a2,nsx,nsy,hx1,hx2,hy1,hy2]=...
update_lag_tail(xu,yu,xv,yv,dS,nls_s,nls_e,xlag,ylag,alpha)
% lagrangian parameters
dxmin=min(diff(xu));
dymin=min(diff(yu));
diagm=sqrt(dxmin^2+dymin^2);
nls1=nls_s-1;
ii=0;
for i=nls_s:nls_e
    ii=ii+1;
    ylag(i)=ylag(nls1)+ii*diagm*sin(alpha);
    xlag(i)=xlag(nls1)+ii*diagm*cos(alpha);
end
[epx,bx,ix1,jy1,a1,nsx,hx1,hy1]=precompute_eps_norkpm(xu(2:length(xu)-1),yu(2:length(yu)-1),xlag,ylag,dS);
[epy,by,ix2,jy2,a2,nsy,hx2,hy2]=precompute_eps_norkpm(xv(2:length(xv)-1),yv(2:length(yv)-1),xlag,ylag,dS);
