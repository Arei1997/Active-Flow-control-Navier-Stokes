function [Uf,Vf]=sumforce(phi,u,v,hx,hy,iwhat)
[ncx,ncy]=size(phi);
nxu=ncx+1; nyu=ncy+2;
nxv=ncx+2; nyv=ncy+1;
coefx=1/hx;
coefy=1/hy;
if iwhat==1
   Uf=zeros(nxu-2,nyu-2); Vf=zeros(nxu-2,nyu-2);
   Uf=u-coefx*(phi(2:ncx,:)-phi(1:ncx-1,:));
   Vf=coefx*(phi(2:ncx,:)-phi(1:ncx-1,:));
else
   Uf=zeros(nxv-2,nyv-2); Vf=zeros(nxv-2,nyv-2);
   Uf=v-coefy*(phi(:,2:ncy)-phi(:,1:ncy-1));
   Vf=coefy*(phi(:,2:ncy)-phi(:,1:ncy-1));
end
