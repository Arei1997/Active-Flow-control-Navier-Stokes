function [Ub]=setBCout(u,ux,uslice,dt,dx,Qin,it,icomp)
[nx,ny]=size(u); [nxu,nyu]=size(ux);
if it > 1
   for j=1:ny
       den=(uslice(2,j)-uslice(1,j));
       if abs(den) > eps
          ratio=(ux(nx-1,j)-uslice(2,j))/den; 
          Uc=-dx/dt*ratio;
          if Uc <= 0
             Uc=0;
          elseif Uc>= dx/dt
             Uc=dx/dt;
          end
       else
          Uc=sum(ux(nxu-1,2:nyu-1))/(nyu-2);
       end
       if icomp==1 & mod(it,10)==0
       disp('Uc')
       disp(Uc)
       end
       Ub(j)=u(nx,j)-dt/dx*Uc*(u(nx,j)-u(nx-1,j));
   end
else
   Uc=sum(ux(nxu-1,2:nyu-1))/(nyu-2);
   for j=1:ny
       Ub(j)=u(nx,j)-dt/dx*Uc*(u(nx,j)-u(nx-1,j));
   end
end
%if icomp==1
%   Qout=sum(Ub(2:ny-1));
%   Ub=Ub*Qin/Qout;
%else
   if icomp==2
   Qout=sum(Ub(2:ny-1));
   Ub=Ub-Qout/(ny-2);
   Ub=zeros(ny,1);
end
