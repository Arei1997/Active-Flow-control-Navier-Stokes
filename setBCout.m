function [Ub]=setBCout(u,v,dt,dx,Qin,icomp)
%function [Ub]=setBCout(u,dt,dx,Qin,icomp) %prev. call
[nx,ny]=size(u); 
Uc=sum(u(nx-1,2:ny-1))/(ny-2);
if icomp==1
   for j=1:ny
       Ub(j)=u(nx,j)-dt/dx*Uc*(u(nx,j)-u(nx-1,j));
   end
  Qout=sum(Ub(2:ny-1));
  Ub=Ub*Qin/Qout;
else
%   Ub=zeros(ny,1);
   [nxv,nyv]=size(v);
   for j=1:nyv
       Ub(j)=v(nxv-1,j)-dt/dx*Uc*(v(nxv-1,j)-v(nxv-2,j));
   end

end
