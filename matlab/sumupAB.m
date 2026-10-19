function [rhs]=sumupAB(u,lapla,NL,laplao,NLo,Gp,Deltat,Re,AB1,AB2)
[nx,ny]=size(u);
rhs=zeros(nx-2,ny-2); 
rhs=-u(2:nx-1,2:ny-1)/Deltat+AB1*(NL-lapla/Re)+AB2*(laplao/Re-NLo)+Gp;
%rhs=-u(2:nx-1,2:ny-1)/Deltat+AB1*NL+5*lapla/Re/2+AB2*(laplao/Re-NLo)+Gp;

