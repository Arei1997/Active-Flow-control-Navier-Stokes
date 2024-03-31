function div=divergence(U,V,hx,hy,ncx,ncy)
div=zeros(ncx,ncy);
div=U(2:ncx+1,2:ncy+1)-U(1:ncx,2:ncy+1);
div=div/hx+(V(2:ncx+1,2:ncy+1)-V(2:ncx+1,1:ncy))/hy;
