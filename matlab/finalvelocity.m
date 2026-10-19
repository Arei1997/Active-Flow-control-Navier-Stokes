function [Uf,Vf,Po]=finalvelocity(phi,u,v,deltat,hx,hy)
[ncx,ncy]=size(phi);
nxu=ncx+1; nyu=ncy+2;
nxv=ncx+2; nyv=ncy+1;
Po=zeros(ncx+2,ncy+2);
Uf=zeros(nxu,nyu); Vf=zeros(nxv,nyv);
for j=1:ncy
    for i=1:ncx
        Po(i+1,j+1)=phi(i,j);
    end
end
for j=2:ncy+1
    Po(1,j)=Po(2,j);
    Po(ncx+2,j)=Po(ncx+1,j);
end
for i=2:ncx+1
    Po(i,1)=Po(i,2);
    Po(i,ncy+2)=Po(i,ncy+1);
end
Po(1,1)=phi(1,1);
Po(1,ncy+2)=phi(1,ncy);
Po(ncx+2,ncy+2)=phi(ncx,ncy);
Po(ncx+2,1)=phi(ncx,1);
for j=1:nyu
    for i=1:nxu
        Uf(i,j)=-deltat*(Po(i+1,j)-Po(i,j))/hx+u(i,j);
    end
end
for j=1:nyv
    for i=1:nxv
        Vf(i,j)=-deltat*(Po(i,j+1)-Po(i,j))/hy+v(i,j);
    end
end
