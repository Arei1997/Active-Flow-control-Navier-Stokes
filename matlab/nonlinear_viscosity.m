function [ustar,vstar]=nonlinear_viscosity(u,v,lx,ly,Deltat,ncx,ncy,Re)
hx=lx/ncx;
hy=ly/ncy;
[nxu,nyu]=size(u);
[nxv,nyv]=size(v);
% A1 is the non linear term for the x-mom equation
for j=2:nyu-1
    for i=2:nxu-1
        A1(i-1,j-1)=(u(i+1,j)^2-u(i-1,j)^2)/(2*hx);
        uavetop=0.5*(u(i,j+1)+u(i,j));
        uavebot=0.5*(u(i,j)+u(i,j-1));
    end 
end
for j=2:nyv-1
    for i=2:nxv-1
        vavetop=0.5*(v(i+1,j)+v(i,j));
        vavebot=0.5*(v(i+1,j-1)+v(i,j-1));
        %A1(i-1,j-1)=A1(i-1,j-1)+(vavetop*uavetop-vavebot*uavebot)/hy;
    end
end
A1=A1+(vavetop*uavetop-vavebot*uavebot)/hy;

% A2 is the non linear term for the y-mom equation
for j=2:nyv-1
    for i=2:nxv-1
        A2(i-1,j-1)=(v(i,j+1)^2-v(i,j-1)^2)/(2*hy);
        vaveleft=0.5*(v(i,j)+v(i-1,j));
        vaveright=0.5*(v(i+1,j)+v(i,j));
    end 
end
for j=2:nyu-1
    for i=2:nxu-1
        uaveleft=0.5*(u(i-1,j+1)+u(i-1,j));
        uaveright=0.5*(u(i,j+1)+u(i,j));
        %A2(i-1,j-1)=A2(i-1,j-1)+(uaveright*vaveright-uaveleft*vaveleft)/hx;
    end
end
A2=A2+(uaveright*vaveright-uaveleft*vaveleft)/hx;


for j=1:nyu-2
    for i=1:nxu-2
        laplau(i+1,j+1)=(u(i+1,j)+u(i+1,j+2)-2*u(i+1,j+1))/hy^2+(u(i,j+1)-2*u(i+1,j+1)+u(i+2,j+1))/hx^2;
        ustar(i,j)=u(i+1,j+1)+Deltat*(-A1(i,j)+1/Re*laplau(i+1,j+1));
    end  
end
% update ustar and vstar
% ustar
%for j=2:nyu-1
%    for i=2:nxu-1
%        ustar(i-1,j-1)=u(i,j);
%    end
%end
%for j=1:nyu-2
%    for i=1:nxu-2
%        ustar(i,j)=u(i+1,j+1)+Deltat*(-A1(i,j)+1/Re*laplau(i,j));
%    end
%end
% vstar
%for j=2:nyv-1
%    for i=2:nxv-1
%        vstar(i-1,j-1)=v(i,j);
%    end
%end
laplav=zeros(nxv,nyv);
for j=1:nyv-2
    for i=1:nxv-2
        laplav(i+1,j+1)=(v(i+1,j)+v(i+1,j+2)-2*v(i+1,j+1))/hy^2+(v(i,j+1)-2*v(i+1,j+1)+v(i+2,j+1))/hx^2;
        vstar(i,j)=v(i+1,j+1)+Deltat*(-A2(i,j)+1/Re*laplav(i+1,j+1));
    end
end

%for j=1:nyv-2
%    for i=1:nxv-2
%        vstar(i,j)=v(i+1,j+1)+Deltat*(-A2(i,j)+1/Re*laplav(i,j));
%    end
%end
