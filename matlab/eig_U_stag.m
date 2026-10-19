function [P,L]=eig_U_stag(n,m,hx,hy,alpha)
A=zeros(n,n);
scal=-2/hx^2-2/hy^2-alpha;
A=A+scal*eye(n,n)+diag(diag(eye(n-1,n-1)),1)/hx^2+diag(diag(eye(n-1,n-1)),-1)/hx^2;
[P,D]=eig(A);
L=diag(D);
