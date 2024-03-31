function [u]=trispecial(diago,b,ifix)
n=length(diago);
% dont fix value ifix=0
alpha=zeros(n,1);
beta=zeros(n,1);
u=zeros(n,1);
x=zeros(n,1);
if ifix== 0
   alpha(1)=diago(1);
   beta(1)=0;
   for i=2:n
       beta(i)=1/alpha(i-1);
       alpha(i)=diago(i)-beta(i);
   end
   x(1)=b(1);
   for i=2:n
       x(i)=b(i)-beta(i)*x(i-1);
   end
   u(n)=x(n)/alpha(n);
   for i=n-1:-1:1
       u(i)=(x(i)-u(i+1))/alpha(i); 
   end
else
   alpha(1)=diago(1);
   beta(1)=0;
   alpha(2)=diago(2);
   beta(2)=1/diago(1);
   for i=3:n
       beta(i)=1/alpha(i-1);
       alpha(i)=diago(i)-beta(i);
   end
   x(1)=0;
   for i=2:n
       x(i)=b(i)-beta(i)*x(i-1);
   end
   u(n)=x(n)/alpha(n);
   for i=n-1:-1:2
       u(i)=(x(i)-u(i+1))/alpha(i);
   end
   u(1)=x(1);
end

