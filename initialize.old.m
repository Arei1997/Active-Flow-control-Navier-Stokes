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
for i=1:ncx+1
    xu(i)=(i-1)*hx;
end
for i=1:ncx+2
    xv(i)=-hx/2+(i-1)*hx;
end
xp=xv(2:ncx+1);
for j=1:ncy+1
    yv(j)=(j-1)*hy;
end
for j=1:ncy+2
    yu(j)=-hy/2+(j-1)*hy;
end
yp=yu(2:ncy+1);
u=zeros(ncx+1,ncy+2);
v=zeros(ncx+2,ncy+1);
%for j=1:ncy+2
%    y=yu(j);
%    for i=1:ncx+1
%        x=xu(i);
%        u(i,j)=sin(pi/2*y)*sin(pi*x);
%    end
%end
%for j=1:ncy+1
%    y=yu(j);
%    for i=1:ncx+2
%        x=xu(i);
%        v(i,j)=sin(pi*y)*sin(pi*x);
%    end
%end
p=zeros(ncx,ncy);
%define the initial conditions on the edges 
%%%%% left edge: i=1, forall j
% for u component
for j=1:ncy+2
    u(1,j)=Ub(4);
end
% for v component
for j=2:ncy
    v(1,j)=2*Vb(4)-v(2,j);
end
%%%%% right edge: i=i_max and forall j
% for u component
for j=1:ncy+2
    u(ncx+1,j)=Ub(2);                            %here
end
% for v component
for j=2:ncy
    v(ncx+2,j)=2*Vb(2)-v(ncx+1,j);
end
%%%%% bottom edge: j=1 and forall i
for i=1:ncx
   % u(i,1)=u(i,2);  
     u(i,1)=Ub(1);                          %here
end
for i=1:ncx+2
    v(i,1)=Vb(1);
end
%%%%% top edge: j=j_max and forall i
for i=1:ncx
   %u(i,ncy+2)=u(i,ncy+1);                          %here
    u(i,ncy+2)=Ub(3);
end
for i=1:ncx+2
    v(i,ncy+1)=Vb(3);
end
xpsi=xu(2:ncx); 
ypsi=yu(2:ncy)-(ly/ncy)/2;
