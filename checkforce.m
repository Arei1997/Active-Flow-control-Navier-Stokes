function dummy=checkforce(phi,phi0,r1,r2,Re,hx,hy)
[ncx,ncy]=size(phi);
nxu=ncx+1; nyu=ncy+2;
nxv=ncx+2; nyv=ncy+1;
Uf=zeros(nxu,nyu); Vf=zeros(nxv,nyv);
coefx=1/hx;
coefy=1/hy;
% update u
Uf(2:nxu-1,2:nyu-1)=r1(2:nxu-1,2:nyu-1)-coefx*(phi(2:ncx,:)-phi(1:ncx-1,:));
Uf(2:nxu-1,1)=r1(2:nxu-1,1)-coefx*(phi(2:ncx,1)-phi(1:ncx-1,1));
Uf(2:nxu-1,nyu)=r1(2:nxu-1,nyu)-coefx*(phi(2:ncx,ncy)-phi(1:ncx-1,ncy));
Uf(1,:)=r1(1,:);
Uf(nxu,:)=r1(nxu,:);
% update v
Vf(2:nxv-1,2:nyv-1)=r2(2:nxv-1,2:nyv-1)-coefy*(phi(:,2:ncy)-phi(:,1:ncy-1));
Vf(1,2:nyv-1)=r2(1,2:nyv-1)-coefy*(phi(1,2:ncy)-phi(1,1:ncy-1));
Vf(nxv,2:nyv-1)=r2(nxv,2:nyv-1)-coefy*(phi(ncx,2:ncy)-phi(ncx,1:ncy-1));
Vf(:,1)=r2(:,1);
Vf(:,nyv)=r2(:,nyv);
