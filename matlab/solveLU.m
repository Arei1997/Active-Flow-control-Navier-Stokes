function [u]=solveLU(alfa,beta,c,f)
% to be used after LU factorization
n=length(f);
y=zeros(n,1);
y(1)=f(1);
for i=2:n
    y(i)=f(i)-beta(i)*y(i-1);
end
u(n)=y(n)/alfa(n);
for i=n-1:-1:1
    u(i)=(y(i)-c(i)*u(i+1))/alfa(i);
end
