function [xu,xv,xp,yu,yv,yp,u,v,p]=initialize(ncx,ncy,lx,ly,Ub,Vb);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%                                                                               %
%               Ub(3) and Vb(3)                                                 %
%           ¡---------------------¡                                             %
%   Ub(4)   ¡                     ¡ Ub(2) and Vb(2)                             % 
%   Vb(4)   ¡                     ¡                                             %
%           -----------------------                                             %
%               Ub(1) and Vb(1)                                                 %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
hx=lx/ncx;
hy=ly/ncy;
nxu=ncx+1;
nyu=ncy+2;
nxv=ncx+2;
nyv=ncy+1;
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
xp=xv(2:ncx+1);
yp=yu(2:ncy+1);
u=ones(nxu,nyu);
%u=zeros(nxu,nyu);
v=zeros(nxv,nyv);
p=zeros(ncx,ncy);
%define the initial conditions on the edges 
%%%%% left edge
% for u component
for j=1:nyu
    u(1,j)=Ub(4,j);
    u(nxu,j)=Ub(4,j);
end
% for v component
for j=1:nyv
    v(1,j)=2*Vb(4,j)-v(2,j);
end
%%%%% right edge
for j=1:nyv
    v(nxv,j)=2*Vb(2,j)-v(nxv-1,j);
end
% bottom edge
for i=1:nxv
    v(i,1)=Vb(1,i);
end
%%%%% top edge
for i=1:nxv
    v(i,nyv)=Vb(3,i);
end
xpsi=xu(2:ncx); 
ypsi=yu(2:ncy)-(ly/ncy)/2;
