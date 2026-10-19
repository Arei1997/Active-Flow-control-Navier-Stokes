function [P,P1,L,EE]=eig_P_stag(n,m,hx,hy)
A=zeros(n,n);
scal=-2/hx^2-2/hy^2;
A=A+scal*eye(n,n)+diag(diag(eye(n-1,n-1)),1)/hx^2+diag(diag(eye(n-1,n-1)),-1)/hx^2;
A(1,1)=-1/hx^2-2/hy^2;
A(n,n)=-1/hx^2-2/hy^2;
[P,D]=eig(A);
L=diag(D);
a=ones(m,1)*L(m);
a(1)=a(1)+1/hy^2;
a(m)=a(m)+1/hy^2;
B=diag(a)+diag(1/hy^2*ones(1,m-1),1)+diag(1/hy^2*ones(1,m-1),-1);
[P1,E]=eig(B); EE=diag(E);
EE(m)=1;
clear B; clear E;

