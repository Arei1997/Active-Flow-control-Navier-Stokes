function [Uf,Vf,Po]=finalvelocityVK(phi,u,v,Pold,deltat,Re,hx,hy)
[ncx,ncy]=size(phi);
nxu=ncx+1; nyu=ncy+2;
nxv=ncx+2; nyv=ncy+1;
Uf=zeros(nxu,nyu); Vf=zeros(nxv,nyv);
coefx=deltat/hx;
coefy=deltat/hy;
% update u
Uf(2:nxu-1,2:nyu-1)=u(2:nxu-1,2:nyu-1)-coefx*(phi(2:ncx,:)-phi(1:ncx-1,:));
Uf(2:nxu-1,1)=u(2:nxu-1,1)-coefx*(phi(2:ncx,1)-phi(1:ncx-1,1));
Uf(2:nxu-1,nyu)=u(2:nxu-1,nyu)-coefx*(phi(2:ncx,ncy)-phi(1:ncx-1,ncy));
Uf(1,:)=u(1,:);
Uf(nxu,:)=u(nxu,:);
% update v
Vf(2:nxv-1,2:nyv-1)=v(2:nxv-1,2:nyv-1)-coefy*(phi(:,2:ncy)-phi(:,1:ncy-1));
Vf(1,2:nyv-1)=v(1,2:nyv-1)-coefy*(phi(1,2:ncy)-phi(1,1:ncy-1));
Vf(nxv,2:nyv-1)=v(nxv,2:nyv-1)-coefy*(phi(ncx,2:ncy)-phi(ncx,1:ncy-1));
Vf(:,1)=v(:,1);
Vf(:,nyv)=v(:,nyv);
% update p
cophi=0.5*deltat/Re;
lphi=laplaphi(phi,hx,hy);
Po=Pold+phi-cophi*lphi;
clear lphi;
