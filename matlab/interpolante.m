function [val]=interpolante(uact,xx,yy,xlag,ylag,dS,b,ix,jy,area,nsup,hx,hy)
nls=length(xlag);
ni=length(xx);
nj=length(yy);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spread the lagrangian forcing on the eulerian field
val=interpol(uact,nls,ni,nj,xlag,ylag,dS,xx,yy,b,ix,jy,area,nsup,hx,hy);
