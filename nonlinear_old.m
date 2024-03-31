function [ustar,vstar,uaverl,check]=nonlinear(u,v,lx,ly,deltat,ncx,ncy,beta)
hx=lx/ncx; hx2=hx^2; twohx=2*hx; 
hy=ly/ncy; hy2=hy^2; twohy=2*hy;
fxnlx=0.25/hy;
fynly=0.25/hx;
[nxu,nyu]=size(u);
[nxv,nyv]=size(v);
A1=zeros(nxu-2,nyu-2);
A2=zeros(nxv-2,nyv-2);
ustar=zeros(nxu-2,nyu-2); vstar=zeros(nxv-2,nyv-2);
% A1 is the non linear term for the x-mom equation
for j=2:nyu-1
    for i=2:nxu-1
        A1(i-1,j-1)=(u(i+1,j)^2-u(i-1,j)^2)/twohx;
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
A1=A1+fxnlx*(vavetop*uavetop-vavebot*uavebot);

% A2 is the non linear term for the y-mom equation
for j=2:nyv-1
    for i=2:nxv-1
        A2(i-1,j-1)=(v(i,j+1)^2-v(i,j-1)^2)/twohy;
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
A2=A2+fynly*(uaveright*vaveright-uaveleft*vaveleft);
% update ustar and vstar
% ustar
%for j=2:nyu-1
%    for i=2:nxu-1
%        ustar(i-1,j-1)=u(i,j);
%    end
%end
for j=1:nyu-2
    for i=1:nxu-2
        ustar(i,j)=u(i+1,j+1)-deltat*A1(i,j);
    end
end
% vstar
%for j=2:nyv-1
%    for i=2:nxv-1
%        vstar(i-1,j-1)=v(i,j);
%    end
%end
for j=1:nyv-2
    for i=1:nxv-2
        vstar(i,j)=v(i+1,j+1)-deltat*A2(i,j);
    end
end
%calculate uaverage 
for j=1:nyu
    uaverage=sum(u(nxu-1,j));
    uaverage=uaverage/nyu;
    uaverl(nxu,j)=u(nxu,j)-uaverage*beta*(u(nxu,j)-u(nxu-1,j));
end
check=sum(uaverl(nxu,j));

