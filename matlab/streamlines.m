function [hlines]=streamlines(u,v,xp,yp)
[nxu,nyu]=size(u);
[nxv,nyv]=size(v);
for j=2:nyu-1
    for i=1:nxu-1
        U(i,j-1)=(u(i+1,j+1)+u(i,j+1))/2;
    end
end

for j=1:nyv-1
    for i=2:nxv-1
        V(i-1,j)=(v(i+1,j)+v(i+1,j+1))/2;
    end
end

[X,Y]=meshgrid(xp,yp);
[sx,sy]=meshgrid(2,-5:5);
hlines=streamline(X,Y,U,V,sx,sy);
set(hlines,'Color','red')
view(2)

