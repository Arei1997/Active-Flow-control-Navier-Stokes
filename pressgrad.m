function [dp]=pressgrad(phi,hx,hy,icomp)
[ncx,ncy]=size(phi);
if icomp==1
   dp=zeros(ncx-1,ncy);
   dp=phi(2:ncx,:)-phi(1:ncx-1,:);
   dp=dp/hx;
elseif icomp==2
   dp=zeros(ncx,ncy-1);
   dp=phi(:,2:ncy)-phi(:,1:ncy-1);
   dp=dp/hy;
end
