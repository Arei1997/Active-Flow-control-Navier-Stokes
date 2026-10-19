function [P,P1,L]=preassemble_stag_P(n,m,hx,hy)

%A=zeros(n,n);
%for i=2:n
%    scal=-(2+2*hy^2/hx^2);
%    A(i,i)=scal;
%end
%A(1,1)=-(hy^2/hx^2+2);
%A(n,n)=-(hy^2/hx^2+2);
%A=A+hy^2*diag(diag(eye(n-1,n-1)),1)/hx^2+hy^2*diag(diag(eye(n-1,n-1)),-1)/hx^2;
%[P,D]=eig(A);
%L=diag(D);
%P1=inv(P);


A=zeros(n,n);
A(1,1)=-(hy^2/hx^2+2);
for i=2:n
    scal=-(2+2*hy^2/hx^2);
    A(i,i)=scal;
end
A(n,n)=-(hy^2/hx^2+2);
A=A+hy^2*diag(diag(eye(n-1,n-1)),1)/hx^2+hy^2*diag(diag(eye(n-1,n-1)),-1)/hx^2;
[P,D]=eig(A);
L=diag(D);
P1=inv(P);
