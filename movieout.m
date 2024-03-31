function [ierr]=movieout(ii,xc,yc,rc,xp,yp,xlag,ylag,u,v)
hx=xp(2)-xp(1); hy=yp(2)-yp(1);
ic=ceil(xc/hx)+1; jc=ceil(yc/hy)+1;
ir1=ceil((xc-rc)/hx); ir2=ceil((xc+rc)/hx);
idi=ir2-ir1;
istart=ic-floor(0.75*idi); iend=ic+floor(6.75*idi);
jstart=jc-floor(1.5*idi); jend=jc+floor(1.5*idi);
a6=axes;
[X,Y]=meshgrid(xp(istart:iend),yp(jstart:jend));
[omega,psi,xo,yo]=vort(u,v,hx,hy);
minval=min(min(omega)); maxval=max(max(omega));
lung=0.8*(maxval-minval); inc=lung/400;
conts=0.8*minval:4*inc:0.8*maxval;
%conts=linspace(minval,maxval,100);
contour(X,Y,(omega(istart:iend,jstart:jend))',conts);
set(a6,'DataAspectRatio',[1 1 1],'XLim',[xp(istart) xp(iend)],'YLim',[yp(jstart) yp(jend)]);
axis equal
title(sprintf(' vorticity  (%g...%g); ',...
                minval,maxval));
colorbar;
hold on; plot(xlag,ylag,'k.');
drawnow;
name_of_file=sprintf('RES/%d.png',ii);
print('-dpng',name_of_file);
hold off;
ierr=1;
