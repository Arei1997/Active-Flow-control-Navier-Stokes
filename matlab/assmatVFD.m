function [PX,PY,EX,EY]=assmatVFD(nxv,nyv,hx,hy)
% use D on all
nx=nxv-2; ny=nyv-2;
A=zeros(nx,nx);
scal=-2/hx^2;
A=A+scal*eye(nx,nx)+diag(diag(eye(nx-1,nx-1)),1)/hx^2+diag(diag(eye(nx-1,nx-1)),-1)/hx^2;
A(1,1)=-3/hx^2;
A(nx,nx)=-3/hx^2;
[PX,D]=eig(A);
EX=diag(D);
clear A; clear D;
A=zeros(ny,ny);
scal=-2/hy^2;
A=A+scal*eye(ny,ny)+diag(diag(eye(ny-1,ny-1)),1)/hy^2+diag(diag(eye(ny-1,ny-1)),-1)/hy^2;
[PY,D]=eig(A);
EY=diag(D);
%EMIX=zeros(nx+1,ny);
%EMIX=ones(nx+1,1)*EY'+EX*ones(ny,1)';
clear A; clear D; 
