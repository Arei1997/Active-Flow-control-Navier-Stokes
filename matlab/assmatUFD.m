function [PX,PY,EX,EY]=assmatUFD(nx,ny,hx,hy)
% OJO it is solved also on right bc
% use D on west; N N N on south north east
A=zeros(nx-2,nx-2);
scal=-2/hx^2;
A=A+scal*eye(nx-2,nx-2)+diag(diag(eye(nx-3,nx-3)),1)/hx^2+...
  diag(diag(eye(nx-3,nx-3)),-1)/hx^2;
[PX,D]=eig(A);
EX=diag(D);
clear A; clear D;
A=zeros(ny-2,ny-2);
scal=-2/hy^2;
A=A+scal*eye(ny-2,ny-2)+diag(diag(eye(ny-3,ny-3)),1)/hy^2+...
  diag(diag(eye(ny-3,ny-3)),-1)/hy^2;
A(1,1)=-1/hy^2; 
A(ny-2,ny-2)=-1/hy^2;
[PY,D]=eig(A);
EY=diag(D);
clear A; clear D; 
