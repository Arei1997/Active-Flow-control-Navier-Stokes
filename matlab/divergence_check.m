function divmax=divergence_check(U,V,hx,hy,ncx,ncy)
for j=1:ncy
    for i=1:ncx
        div(i,j)=(U(i+1,j+1)-U(i,j+1))/hx+(V(i+1,j+1)-V(i+1,j))/hy;
    end
end
divmax=max(max(abs(div)));
clear div;
