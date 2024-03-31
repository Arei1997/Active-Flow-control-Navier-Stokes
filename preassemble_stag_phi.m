function [P,P1,L]=preassemble_stag_phi(n,m,hx,hy,alpha)
A=zeros(n,n);
scal=-(2*+2*hy^2/hx^2+alpha*hy^2);
A=A+scal*eye(n,n)+hy^2*diag(diag(eye(n-1,n-1)),1)/hx^2+hy^2*diag(diag(eye(n-1,n-1)),-1)/hx^2;
[P,D]=eig(A);
L=diag(D);
P1=inv(P);
%
