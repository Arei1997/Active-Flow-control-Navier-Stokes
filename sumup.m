function [rhs]=sumup(u,lapla,NL,Deltat,Re,icomp)
[nx,ny]=size(u);
% A1 is the non linear term for the x-mom equation
if icomp==1
rhs=zeros(nx-2,ny-2); 
    for j=2:ny-1
        for i=2:nx-1
            rhs(i-1,j-1)=-u(i,j)/Deltat-lapla(i-1,j-1)/Re+NL(i-1,j-1);
        end
     end
elseif icomp==2
     rhs=zeros(nx-2,ny-2); 
     for j=2:ny-1
         for i=2:nx-1
             rhs(i-1,j-1)=-u(i,j)/Deltat-lapla(i-1,j-1)/Re+NL(i-1,j-1);
         end
     end
end

