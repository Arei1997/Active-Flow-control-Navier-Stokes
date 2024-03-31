function [uoo,voo]=setbc(Uoo,Voo,Vb,Ub)
[n,m]=size(Uoo);
nxu=n+2;
nyu=m+2;
[n,m]=size(Voo);
nxv=n+2;
nyv=m+2;
ncx=nxu-1; ncy=nyv-1;
uoo=zeros(nxu,nyu);
voo=zeros(nxv,nyv);
for i=2:nxu-1
    for j=2:nyu-1
        uoo(i,j)=Uoo(i-1,j-1);
    end
end
for i=2:nxv-1
    for j=2:nyv-1
        voo(i,j)=Voo(i-1,j-1);
    end
end
%%%%% left edge: i=1, forall j
% for u component
for j=1:ncy+2
    uoo(1,j)=Ub(4);
end
% for v component
for j=2:ncy
    voo(1,j)=2*Vb(4)-voo(2,j);
end
%%%%% right edge: i=i_max and forall j
% for u component
for j=1:ncy+2
    %uaverage=sum(uoo(ncx,j))/(ncy+2);
    uoo(ncx+1,j)=Ub(2);
    %uoo(ncx+1,j)=uoo(ncx,j);                            %here
end
% for v component
for j=2:ncy
    voo(ncx+2,j)=2*Vb(2)-voo(ncx+1,j);
end
%%%%% bottom edge: j=1 and forall i
for i=2:ncx
    %uoo(i,1)=2*Ub(1)-uoo(i,2);
    uoo(i,1)=uoo(i,2);                                %here
end
for i=1:ncx+2
    voo(i,1)=Vb(1);
end
%%%%% top edge: j=j_max and forall i
for i=2:ncx
    %uoo(i,ncy+2)=2*Ub(3)-uoo(i,ncy+1);
   uoo(i,ncy+2)=uoo(i,ncy+1);                          %here
end
for i=1:ncx+2
    voo(i,ncy+1)=Vb(3);
end


