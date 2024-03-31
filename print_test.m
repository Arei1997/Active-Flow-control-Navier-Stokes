function print_test(T,X,Y,xlag,ylag,nj,ni)

t=reshape(T,nj-2,ni-2);
pcolor(X,Y,t); axis equal; hold on;title('campo corregido');plot(xlag,ylag,'o-r');
%plot(xe,ye,'.r');plot(xi,yi,'.b');
%figure(10);contourf(X,Y,t,v); shading faceted; axis equal; hold on; plot(xlag,ylag,'o-r');

caxis([-1e4 1e4]);colormap(jet2);

% pcolor(X,Y,t); axis equal; hold on; plot(xlag,ylag,'w*');pause(0.01);hold off;shading faceted
drawnow;pause(0.1); hold off;shading interp;