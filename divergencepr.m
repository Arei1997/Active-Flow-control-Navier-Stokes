function div=divergencepr(U,V,hx,hy,ncx,ncy,deltat)
div=zeros(ncx,ncy);
for j=1:ncy
    for i=1:ncx
        div(i,j)=1/deltat*(hy^2*(U(i+1,j+1)-U(i,j+1))/hx+hy*(V(i+1,j+1)-V(i+1,j)));
    end
end
