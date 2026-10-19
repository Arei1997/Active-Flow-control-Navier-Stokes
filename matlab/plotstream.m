clear all
load uu;
load vv;
load xp;
load yp;
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
xmax=max(xp); xmin=min(xp); xpsize=length(xp);
ymax=max(yp); ymin=min(yp); ypsize=length(yp);
disp('domain in x: xmax,xmin,nxp'); disp(xmax); disp(xmin); disp(xpsize);
disp('domain in y: ymax,ymin,nyp'); disp(ymax); disp(ymin); disp(ypsize);
itimes=0;
iok=0;
while iok==0
   itimes=itimes+1;
   is=input('enter first index in x= ');
   ie=input('enter last index in x= ');
   js=input('enter first index in y= ');
   je=input('enter last index in y= ');
   ist=input('enter stride in x= ');
   jst=input('enter stride in y= ');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if itimes >1
      clear uplot; clear vplot; clear X; clear Y;
      close(f1);
   end
   uplot=u(is:ist:ie,js:jst:je); uplot=uplot';
   vplot=v(is:ist:ie,js:jst:je); vplot=vplot';
   [X,Y]=meshgrid(xp(is:ist:ie),yp(js:jst:je));
   vecsc=input('scale for vector plot= ');
   f1=figure;
   quiver(X,Y,uplot,vplot,vecsc);
   axis equal;
   iok=input('if plot needs to be redone answer 0 (any other number if ok) ');
end
%plot circle
xc=input('center of circle x coordinate= '); 
yc=input('center of circle y coordinate= '); 
diam=input('diameter of circle= '); rcirc=diam/2;
f2=figure;
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
contourf(XM,YM,circle');
axis equal
% plot streamlines
load xu; load yu;
[xlag,ylag,theta,dS]=initialize_lag(xu,yu,rcirc,xc,yc,1,0);
%plot(xlag,ylag,'k*');
istream=floor(0.5*(ie-is)/ist);
xpts = linspace(xp(is),xp(ie),istream);
jstream=floor(0.5*(je-js)/jst);
ypts = linspace(yp(js),yp(je),jstream);
[sx,sy] = meshgrid(xpts,ypts);
stepsize=xp(2)-xp(1);
novert=ie-is;
opstream=[stepsize, novert];
streamline(xp(is:ist:ie),yp(js:jst:je),uplot,vplot,sx,sy,opstream)
