function [lapla]=laplacian(un,hx,hy,Ub,icomp)
% icomp=1 U
% icomp=2 V
[nx1,ny1]=size(un); nx=nx1-2; ny=ny1-2;
u=zeros(nx,ny); u=un(2:nx1-1,2:ny1-1);
lapla=zeros(nx,ny);
for j=2:ny-1
    for i=2:nx-1
        lapla(i,j)=(u(i,j-1)+u(i,j+1)-2*u(i,j))/hy^2+(u(i-1,j)-2*u(i,j)+u(i+1,j))/hx^2;
    end
end
if icomp==1
   for j=2:ny-1
       lapla(1,j)=(Ub(4)-2*u(1,j)+u(2,j))/hx^2+(u(1,j-1)-2*u(1,j)+u(1,j+1))/hy^2; 
       lapla(nx,j)=(u(nx-1,j)-1*u(nx,j))/hx^2+(u(nx,j-1)-2*u(nx,j)+u(nx,j+1))/hy^2; %here 1
   end
   for i=2:nx-1
       % valueb=(2*Ub(1)-u(i,1));
        lapla(i,1)=(-1*u(i,1)+u(i,2))/hy^2+(u(i-1,1)-2*u(i,1)+u(i+1,1))/hx^2;  %here 
       % valuet=2*Ub(3)-u(i,ny);
        lapla(i,ny)=(u(i,ny-1)-1*u(i,ny))/hy^2+(u(i-1,ny)-2*u(i,ny)+u(i+1,ny))/hx^2; %here
   end
   %valueb=2*Ub(1)-u(1,1);
   valueb=u(1,1);
   lapla(1,1)=(valueb-2*u(1,1)+u(1,2))/hy^2+(Ub(4)-2*u(1,1)+u(2,1))/hx^2; %here
   %valueb=2*Ub(1)-u(nx,1);
   valueb=u(nx,1);
   lapla(nx,1)=(valueb-2*u(nx,1)+u(nx,2))/hy^2+(u(nx-1,1)-1*u(nx,1))/hx^2; %here 2
   %valuet=2*Ub(3)-u(1,ny);
   valuet=u(1,ny);
   lapla(1,ny)=(u(1,ny-1)-2*u(1,ny)+valuet)/hy^2+(Ub(4)-2*u(1,ny)+u(2,ny))/hx^2; %here
   %valuet=2*Ub(3)-u(nx,ny);
   valuet=u(nx,ny);
   lapla(nx,ny)=(u(nx,ny-1)-2*u(nx,ny)+valuet)/hy^2+(u(nx-1,ny)-1*u(nx,ny))/hx^2; %here 3
else
   for j=2:ny-1
       valuel=2*Ub(4)-u(1,j);
       lapla(1,j)=(u(1,j-1)-2*u(1,j)+u(1,j+1))/hy^2+(valuel-2*u(1,j)+u(2,j))/hx^2;
       valuer=2*Ub(2)-u(nx,j);
       lapla(nx,j)=(u(nx,j-1)-2*u(nx,j)+u(nx,j+1))/hy^2+(u(nx-1,j)-2*u(nx,j)+valuer)/hx^2;
   end
   for i=2:nx-1
       lapla(i,1)=(Ub(1)-2*u(i,1)+u(i,2))/hy^2+(u(i-1,1)-2*u(i,1)+u(i+1,1))/hx^2;
       lapla(i,ny)=(u(i,ny-1)-2*u(i,ny)+Ub(3))/hy^2+(u(i-1,ny)-2*u(i,ny)+u(i+1,ny))/hx^2;
   end
   valuel=2*Ub(4)-u(1,1);
   lapla(1,1)=(Ub(1)-2*u(1,1)+u(1,2))/hy^2+(valuel-2*u(1,1)+u(2,1))/hx^2;
   valuer=2*Ub(2)-u(nx,1);
   lapla(nx,1)=(Ub(1)-2*u(nx,1)+u(nx,2))/hy^2+(u(nx-1,1)-2*u(nx,1)+valuer)/hx^2;
   valuelu=2*Ub(4)-u(1,ny);
   lapla(1,ny)=(u(1,j-1)-2*u(1,ny)+Ub(3))/hy^2+(valuelu-2*u(1,ny)+u(2,ny))/hx^2;
   valueru=2*Ub(2)-u(nx,ny);
   lapla(nx,ny)=(u(nx,ny-1)-2*u(nx,ny)+Ub(3))/hy^2+(u(nx-1,ny)-2*u(nx,ny)+valueru)/hx^2;
end


