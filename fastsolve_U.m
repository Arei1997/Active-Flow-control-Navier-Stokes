function f=fastsolve_U(P,P1,L,f,Ub,hx,hy,alpha)
[n,m]=size(f);
nyu=m+2; nxu=n+2;
FT=P1*f*hy^2;
UT=zeros(n,m);
%di=zeros(m,1);
for i=1:n
    di(1)=L(i)+1;
    di(2:m-1)=L(i);
    di(m)=L(i)+1;
    xt=zeros(m,1);
    a=zeros(m,1);
    b=zeros(m-1,1);
    c=zeros(m-1,1);
    a(1)=L(i);
    %r(1)=sqrt(L(i));
    for k=1:m-1
        b(k)=1/a(k);
        a(k+1)=L(i)-b(k);
    end
    xt(1)=FT(i,1);
    for k=2:m
        xt(k)=(FT(i,k)-b(k-1)*xt(k-1));
    end
    UT(i,m)=xt(m)/a(m);
    for k=m-1:-1:1
        UT(i,k)=(xt(k)-UT(i,k+1))/a(k);
    end
end
f=P*UT;
