function [NL]=nonlinear(u,v,Ub,Vb,lx,ly,ncx,ncy,icomp)
hx=lx/ncx; twohx=2*hx;
hy=ly/ncy; twohy=2*hy;
fxnlx=0.25/hy;
fynly=0.25/hx;
[nxu,nyu]=size(u);
[nxv,nyv]=size(v);
if icomp==1
   NL=zeros(nxu-2,nyu-2); 
   NL=(v(3:nxu,2:nyu-1)+v(2:nxu-1,2:nyu-1)).*(u(2:nxu-1,3:nyu)+u(2:nxu-1,2:nyu-1));
   NL=NL-(v(3:nxu,1:nyu-2)+v(2:nxu-1,1:nyu-2)).*(u(2:nxu-1,2:nyu-1)+u(2:nxu-1,1:nyu-2));
   NL=fxnlx*NL+(u(3:nxu,2:nyu-1).^2-u(1:nxu-2,2:nyu-1).^2)/twohx;
elseif icomp==2
%%%%%%%%%%%%%%%%%%
   NL=zeros(nxv-2,nyv-2);
   NL=(u(2:nxv-1,3:nyv)+u(2:nxv-1,2:nyv-1)).*(v(2:nxv-1,2:nyv-1)+v(3:nxv,2:nyv-1));
   NL=NL-(u(1:nxv-2,3:nyv)+u(1:nxv-2,2:nyv-1)).*(v(2:nxv-1,2:nyv-1)+v(1:nxv-2,2:nyv-1));
   NL=fynly*NL+(v(2:nxv-1,3:nyv).^2-v(2:nxv-1,1:nyv-2).^2)/twohy;
else
    disp('error in nonlin')
end

