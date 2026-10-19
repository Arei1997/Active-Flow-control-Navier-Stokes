function [u]=specialbid(alpha,beta,f)
% to be used after call to specialLU
% alpha and beta outputs of specialLU
% f rhs of the system; u the solution
n=length(alpha);
c=beta(1);
y=zeros(n,1);
u=zeros(n,1);
y(1)=f(1);
for i=2:n
    y(i)=f(i)-beta(i)*y(i-1);
end
u(n)=y(n)/alpha(n);
for i=n-1:-1:1
    u(i)=(y(i)-c*u(i+1))/alpha(i);
end
