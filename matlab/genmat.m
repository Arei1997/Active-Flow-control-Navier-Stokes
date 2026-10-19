function [DD]=genmat(x,y)
nx=length(x); ny=length(y);
D=zeros(nx,nx);

for i=2:nx-1
    dxm=(x(i)-x(i-1));
    dxp=(x(i+1)-x(i));
    factor=2/(dxp+dxm);
    D(i,i-1)=factor/dxm;
    D(i,i)=-factor*(1/dxm+1/dxp); 
    D(i,i+1)=factor/dxp;
end
Dx=D(2:nx-1,2:nx-1); clear D;
D=zeros(ny,ny);

for i=2:ny-1
    dxm=(y(i)-y(i-1));
    dxp=(y(i+1)-y(i));
    factor=2/(dxp+dxm);
    D(i,i-1)=factor/dxm;
    D(i,i)=-factor*(1/dxm+1/dxp);
    D(i,i+1)=factor/dxp;
end
Dy=D(2:ny-1,2:ny-1); clear D;
DD = -kron(speye(nx-2,nx-2),Dy)-kron(Dx,speye(ny-2,ny-2));

