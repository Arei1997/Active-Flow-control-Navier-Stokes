function [ustar,vstar]=nonlinear_viscosity_new(u,v,lx,ly,Deltat,ncx,ncy,Re)
hx=lx/ncx; hx2=hx^2; twohx=2*hx; 
hy=ly/ncy; hy2=hy^2; twohy=2*hy;
fxnlx=0.25/hy;
fynly=0.25/hx;
[nxu,nyu]=size(u);
[nxv,nyv]=size(v);
% A1 is the non linear term for the x-mom equation
ustar=zeros(nxu-2,nyu-2); vstar=zeros(nxv-2,nyv-2);

for j=2:nyu-1
    for i=2:nxu-1
        nonlin=(u(i+1,j)^2-u(i-1,j)^2)/twohx;
        vavetop=v(i+1,j)+v(i,j);
        vavebot=v(i+1,j-1)+v(i,j-1);
        uavetop=u(i,j+1)+u(i,j);
        uavebot=u(i,j)+u(i,j-1);
        nonlin=nonlin+fxnlx*(vavetop*uavetop-vavebot*uavebot);
        laplau=(u(i,j-1)+u(i,j+1)-2*u(i,j))/hy2+(u(i-1,j)-2*u(i,j)+u(i+1,j))/hx2;
        ustar(i-1,j-1)=u(i,j)+Deltat*(1/Re*laplau-nonlin);
    end
end


%%%%%%%%%%%%%%%%%% 

for j=2:nyv-1
    for i=2:nxv-1
        nonlin=(v(i,j+1)^2-v(i,j-1)^2)/twohy;
        uaveleft=0.5*(u(i-1,j+1)+u(i-1,j));
        uaveright=0.5*(u(i,j+1)+u(i,j));
        vaveleft=0.5*(v(i,j)+v(i-1,j));
        vaveright=0.5*(v(i+1,j)+v(i,j));
        nonlin=nonlin+fynly*(uaveright*vaveright-uaveleft*vaveleft);
        laplav=(v(i,j-1)+v(i,j+1)-2*v(i,j))/hy2+(v(i-1,j)-2*v(i,j)+v(i+1,j))/hx2;
        vstar(i-1,j-1)=v(i,j)+Deltat*(1/Re*laplav-nonlin);
    end
end
