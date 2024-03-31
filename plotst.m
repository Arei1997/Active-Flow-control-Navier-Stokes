clear all
load uu;
load vv;
load xp;
load yp;
hx=xp(2)-xp(1); hy=yp(2)-yp(1);
[omega,psi,xo,yo]=vort(u,v,hx,hy); 
[nx1,ny1]=size(psi);
for i=1:nx1
    psi(i,:)=psi(i,:)+yo(2:length(yo)-1);
end
%plot circle
xc=input('center of circle x coordinate= '); 
yc=input('center of circle y coordinate= '); 
diam=input('diameter of circle= '); rcirc=diam/2;
f2=figure;
%%%% OJO
ie=280; is=40;
je=270; js=150;
xm=linspace(xp(is),xp(ie),10*(ie-is));
ym=linspace(yp(js),yp(je),10*(je-js));
n1=length(xm); n2=length(ym);
circle=zeros(n1,n2);
for i=1:n1
    for j=1:n2
        dist=(xm(i)-xc)^2+(ym(j)-yc)^2;
        if dist <= rcirc^2
             circle(i,j)=1;
        end
    end
end
[XM,YM]=meshgrid(xm,ym);
%contourf(XM,YM,circle');
axis equal
hold on;
% plot streamlines
range=input('range of streamfunction= '); 
nom=input('number of streamlines= '); 
range=range/2;
conts=linspace(yc-range,yc+range,nom);
[X,Y]=meshgrid(xo(2:length(xo)-1),yo(2:length(yo)-1));
contour(X,Y,psi',conts);
