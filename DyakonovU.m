function [us]=DyakonovU(f,un,Ub,Ubs,alx,bex,cx,aly,bey,cy,hx,hy,Re,Deltat)
[m,n]=size(f); g=zeros(m,n); us=zeros(m,n);
nx=m+2; ny=n+2;
alphax=0.5*Deltat/(Re*hx^2);
alphay=0.5*Deltat/(Re*hy^2);
un(1,1)=(hx*un(2,1)+hy*un(1,2))/(hx+hy);
un(nx,1)=(hx*un(nx-1,1)+hy*un(nx,2))/(hx+hy);
un(nx,ny)=(hx*un(nx-1,ny)+hy*un(nx,ny-1))/(hx+hy);
un(1,ny)=(hx*un(2,ny)+hy*un(1,ny-1))/(hx+hy);
g1=zeros(nx,ny-2);
g1=un(:,3:ny)-2*un(:,2:ny-1)+un(:,1:ny-2); g1=alphay*g1+un(:,2:ny-1);
g=g1(3:nx,:)-2*g1(2:nx-1,:)+g1(1:nx-2,:); g=alphax*g+g1(2:nx-1,:);
clear g1;
g=g+Deltat*f;
ubleft=zeros(nx-2,1);
ubright=zeros(nx-2,1);
for j=2:ny-1
    ubleft(j-1)=Ub(4,j)-alphay*(Ub(4,j+1)-2*Ub(4,j)+Ub(4,j-1));
    ubright(j-1)=Ubs(j)-alphay*(Ubs(j+1)-2*Ubs(j)+Ubs(j-1));
end
gg=zeros(nx-2,1);
for j=1:ny-2
    gg=g(:,j);
    gg(1)=gg(1)+alphax*ubleft(j);
    gg(nx-2)=gg(nx-2)+alphax*ubright(j);
    g(:,j)=solveLU(alx,bex,cx,gg);
end
gg=zeros(ny-2,1);
for i=1:nx-2
    gg=g(i,:);
    us(i,:)=solveLU(aly,bey,cy,gg);
end
clear g;
        
