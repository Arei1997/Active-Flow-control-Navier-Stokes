function Po=fastsolver_P(P,P1,L,div)
[n,m]=size(div);
FT=P1*div;
PT=zeros(n,m);
%di=zeros(m,1);
i=n;
di(1)=L(i)-1;
di(2:m-1)=L(i);
di(m)=L(i)-1;
xt=zeros(m,1);
aa=zeros(m,1);
bb=zeros(m-1,1);
cc=zeros(m-1,1);
aa(1)=L(i);
cc(1)=1;
bb(1)=0;
cc(2)=0;
aa(2)=1;
for k=2:m-1
    bb(k)=1/aa(k);
    aa(k+1)=L(i)-cc(k)*bb(k);
    cc(k+1)=1;
end
xt(1)=FT(i,1);
FT(i,2)=0;
xt(2)=FT(i,2);
for k=3:m
    xt(k)=(FT(i,k)-bb(k-1)*xt(k-1));
end
PT(i,m)=xt(m)/aa(m);
for k=m-1:-1:1
    PT(i,k)=(xt(k)-cc(m-1)*PT(i,k+1))/aa(k);
end
for i=1:n-1
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
    VT(i,m)=xt(m)/a(m);
    for k=m-1:-1:1
        VT(i,k)=(xt(k)-VT(i,k+1))/a(k);
    end
end
Po=P*PT;
        
