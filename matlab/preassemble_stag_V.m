function [P,P1,L]=preassemble_stag_V(n,m,hx,hy,alpha)
A=zeros(n,n);
for i=1:n
    scal=-(2+3*hy^2/hx^2+alpha*hy^2);
    A(i,i)=scal;
end
A=A+hy^2*diag(diag(eye(n-1,n-1)),1)/hx^2+hy^2*diag(diag(eye(n-1,n-1)),-1)/hx^2;
A(1,1)=-(2+3*hy^2/hx^2+alpha*hy^2);
A(n,n)=-(2+3*hy^2/hx^2+alpha*hy^2);
[P,D]=eig(A);
L=diag(D);
P1=inv(P);
