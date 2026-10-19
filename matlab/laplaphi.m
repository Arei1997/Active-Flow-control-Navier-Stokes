function [lapla]=laplaphi(uu,hx,hy)
[n1,n2]=size(uu); 
nx=n1+2;
ny=n2+2;
u=zeros(nx,ny);
u(2:nx-1,2:ny-1)=uu;
u(2:nx-1,1)=uu(:,1);
u(2:nx-1,ny)=uu(:,n2);
u(1,2:ny-1)=uu(1,:);
u(nx,2:ny-1)=uu(n1,:);
coef=2/hx^2+2/hy^2;
lapla=zeros(nx-2,ny-2);
lapla=u(3:nx,2:ny-1)+u(1:nx-2,2:ny-1);
lapla=lapla/hx^2+(u(2:nx-1,3:ny)+u(2:nx-1,1:ny-2))/hy^2;
lapla=lapla-coef*u(2:nx-1,2:ny-1);
clear u;
