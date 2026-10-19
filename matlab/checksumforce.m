function [Uf,Vf]=checksumforce(phi,u,v,hx,hy)
[ncx,ncy]=size(phi);
nxu=ncx+1; nyu=ncy+2;
nxv=ncx+2; nyv=ncy+1;
Uf=zeros(nxu,nyu); Vf=zeros(nxv,nyv);
c1=1/hx; c2=1/hy;
% update u
Uf(2:nxu-1,2:nyu-1)=u(2:nxu-1,2:nyu-1)-c1*(phi(2:ncx,:)-phi(1:ncx-1,:));
% update v
Vf(2:nxv-1,2:nyv-1)=v(2:nxv-1,2:nyv-1)-c2*(phi(:,2:ncy)-phi(:,1:ncy-1));
