function [eul_f,val]=immersed_boundary(uact,xx,yy,xlag,ylag,dS,epsilon,b,ix,jy,area,nsup,hx,hy,ns,ne,vel)
nls=length(xlag);
ni=length(xx);
nj=length(yy);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spread the lagrangian forcing on the eulerian field
val=interpol(uact,nls,ni,nj,xlag,ylag,dS,xx,yy,b,ix,jy,area,nsup,hx,hy);
if nargin==17
val(ns:ne)=val(ns:ne)+vel';
end
eul_f=spread(val,epsilon,nls,ni,nj,xlag,ylag,dS,xx,yy,b,ix,jy,area,nsup,hx,hy);
