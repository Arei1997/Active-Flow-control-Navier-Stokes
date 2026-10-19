function [Ub,Vb]=set_velocity_bc(ncx,ncy,lx,ly);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Ub(3)  Vb(3)                            %
%         .---------------------.                       %
%         ¡                     ¡                       %
% Ub(4)   ¡                     ¡ Ub(2)and Vb(2)        %
% Vb(4)   ¡                     ¡                       %
%         -----------------------                       %
%               Ub(1)  Vb(1)                            %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
hx=lx/ncx;
hy=ly/ncy;
nxu=ncx+1;
nyu=ncy+2;
nxv=ncx+2;
nyv=ncy+1;
%
% set bcs for U
     Ub=zeros(4,max(nxu,nyu));
% south and north
     for i=1:nxu
        Ub(1,i)=1;
        Ub(3,i)=1;
     end
% east and west
     for j=1:nyu
         Ub(2,j)= 1;
         Ub(4,j)= 1;
     end
% set bcs for V
     Vb=zeros(4,max(nxv,nyv));
% south and north
     for i=1:nxv
        Vb(1,i)= 0;
        Vb(3,i)= 0;
     end
% east and west
     for j=1:nyv
         Vb(2,j)= 0;
         Vb(4,j)= 0;
     end
