function FT=pressurefix(P,P1,EI,f)
[n,m]=size(f);
FT=zeros(n,m); UT=zeros(n,m);
FT=P1*f;
for i=1:n-1
    diago=EI(i)*ones(m,1);
    diago(1)=diago(1)+1;
    diago(m)=diago(m)+1;
    UT(i,:)=trispecial(diago,FT(i,:),0);
end
diago=EI(n)*ones(m,1);
diago(1)=1;
diago(m)=diago(m)+1;
UT(n,:)=trispecial(diago,FT(n,:),1);
FT=P*UT;
clear UT;


