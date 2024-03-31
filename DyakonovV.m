function [us]=DyakonovV(f,un,Ub,Ubs,alx,bex,cx,aly,bey,cy,hx,hy,Re,Deltat)
%function [us]=DyakonovV(f,un,Ub,alx,bex,cx,aly,bey,cy,hx,hy,Re,Deltat)
[m,n]=size(f); g=zeros(m,n); us=zeros(m,n);
nx=m+2; ny=n+2;
alphax=0.5*Deltat/(Re*hx^2);
alphay=0.5*Deltat/(Re*hy^2);
un(1,1)=(hx*un(2,1)+hy*un(1,2))/(hx+hy);
un(nx,1)=(hx*un(nx-1,1)+hy*un(nx,2))/(hx+hy);
un(nx,ny)=(hx*un(nx-1,ny)+hy*un(nx,ny-1))/(hx+hy);
un(1,ny)=(hx*un(2,ny)+hy*un(1,ny-1))/(hx+hy);
g1=zeros(nx-2,ny);
g1=un(3:nx,:)-2*un(2:nx-1,:)+un(1:nx-2,:); g1=alphax*g1+un(2:nx-1,:);
g=g1(:,3:ny)-2*g1(:,2:ny-1)+g1(:,1:ny-2); g=alphay*g+g1(:,2:ny-1);
clear g1;
g=g+Deltat*f;
ubbot=zeros(nx-2,1);
ubtop=zeros(nx-2,1);
for i=2:nx-1
    %ubbot(i-1)=Ub(1,i)-alphax*(Ub(1,i+1)-2*Ub(1,i)-Ub(1,i-1));
    %ubtop(i-1)=Ub(3,i)-alphax*(Ub(3,i+1)-2*Ub(3,i)-Ub(3,i-1));    
    ubbot(i-1)=Ub(1,i)-alphax*Ub(1,i-1)+2*alphax*Ub(1,i)-alphax*Ub(1,i+1);
    ubtop(i-1)=Ub(3,i)-alphax*Ub(3,i-1)+2*alphax*Ub(3,i)-alphax*Ub(3,i+1);
end
for i=1:nx-2
    gg=zeros(ny-2,1);
    gg=g(i,:);
    gg(1)=gg(1)+alphay*ubbot(i);
    gg(ny-2)=gg(ny-2)+alphay*ubtop(i);
    g(i,:)=solveLU(aly,bey,cy,gg);
end
for j=1:ny-2
    gg=zeros(nx-2,1);
    gg=g(:,j);
    gg(1)=gg(1)+2*alphax*Ub(4,j+1);
%    gg(nx-2)=gg(nx-2)+2*alphax*Ub(2,j+1);
    gg(nx-2)=gg(nx-2)+2*alphax*Ubs(j+1);
    us(:,j)=solveLU(alx,bex,cx,gg);
end
clear g;
        
