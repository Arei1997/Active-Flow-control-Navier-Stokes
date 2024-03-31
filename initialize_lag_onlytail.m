%-------------------------------------------------------------------------
function [xlag,ylag,dS,nls1,nls]=initialize_lag_onlytail(xx,yy,rc,xc,yc,crot,alpha)
% rc is half the length of the beam
% lagrangian parameters
nls1=1;
dxmin=min(diff(xx));
dymin=min(diff(yy));
diagm=sqrt(dxmin^2+dymin^2);
% suppose the tail length is D=2*rc
nls=ceil(2*rc/diagm);
%nls1=nls;
lbeam=2*rc;
nodes=linspace(0,lbeam,nls);
DS_base=nodes(2)-nodes(1);
for i=1:nls
    xlag(i)=xc-crot*cos(alpha)+nodes(i)*cos(alpha);
    ylag(i)=yc-crot*sin(alpha)+nodes(i)*sin(alpha);
    dS(i)=DS_base;
end

