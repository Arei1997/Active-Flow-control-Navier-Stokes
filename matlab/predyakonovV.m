function [alx,bex,cx,aly,bey,cy]=predyakonovV(nx,ny,hx,hy,dt,re)
alphax=0.5*dt/(re*hx^2);
alphay=0.5*dt/(re*hy^2);
alx=zeros(nx-2,1);
bex=zeros(nx-2,1);
cx=zeros(nx-2,1);
ax=zeros(nx-2,1);
ex=zeros(nx-2,1);
for i=1:nx-2
  ax(i)=1+2*alphax;
  ex(i)=-alphax;
  cx(i)=-alphax;
end
ax(1)=1+3*alphax;
ax(nx-2)=1+3*alphax;
[alx,bex]=tridLU(ax,ex,cx);
aly=zeros(ny-2,1);
bey=zeros(ny-2,1);
cy=zeros(ny-2,1);
ay=zeros(ny-2,1);
ey=zeros(ny-2,1);
for i=1:ny-2
  ay(i)=1+2*alphay;
  ey(i)=-alphay;
  cy(i)=-alphay;
end
[aly,bey]=tridLU(ay,ey,cy);
