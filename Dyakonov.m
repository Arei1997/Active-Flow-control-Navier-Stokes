function [g]=Dyakonov(f,un,Ub,alx,bex,cx,aly,bey,cy,hx,hy,Re,Deltat)
[m,n]=size(f); g=zeros(m,n);
nx=m+2; ny=n+2;
alphax=0.5*Deltat/(Re*hx^2);
alphay=0.5*Deltat/(Re*hy^2);
g1=zeros(nx,ny-2);
for i=1:nx
    for j=2:ny-1
        g1(i,j-1)=alphay*(un(i,j+1)-2*un(i,j)+un(i,j-1))+un(i,j);
    end
end
for j=2:ny-1
    for i=2:nx-1
        g(i-1,j-1)=alphax*(g1(i+1,j-1)-2*g1(i,j-1)+g1(i-1,j-1))+g1(i,j-1);
    end
end
clear g1;
g=g+Deltat*f;
ubleft=zeros(nx-2,1);
ubright=zeros(nx-2,1);
for j=2:ny-1
    ubleft(j-1)=Ub(4,j)-alphay*(Ub(4,j+1)-2*Ub(4,j)-Ub(4,j-1));
    ubright(j-1)=Ub(2,j)-alphay*(Ub(2,j+1)-2*Ub(2,j)-Ub(2,j-1));
end
for j=1:ny-2
    gg=zeros(nx-2,1);
    gg=g(:,j);
    gg(1)=gg(1)+alphax*ubleft(j);
    gg(nx-2)=gg(nx-2)+alphax*ubright(j);
    g(:,j)=solveLU(alx,bex,cx,gg);
end
for i=1:nx-2
    gg=zeros(ny-2,1);
    gg=g(i,:);
    g(i,:)=solveLU(aly,bey,cy,gg);
end
        
