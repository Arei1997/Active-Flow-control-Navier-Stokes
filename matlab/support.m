function [ix,jy,area,nsup,hx,hy]=support(xt,yt,x,y,ni,nj)

dmin2=1e10;
for i=1:ni
    for j=1:nj
         d2=(x(i)-xt)^2+(y(j)-yt)^2;
           if d2 < dmin2
            dmin2=d2;
            ii=i;
            jj=j;
        end
    end
end


dx=abs(xt-x(ii)); dy=abs(yt-y(jj));
hxp=max(abs(x(ii+1)-x(ii)),abs(x(ii)-x(ii-1)));
hxm=min(abs(x(ii+1)-x(ii)),abs(x(ii)-x(ii-1)));
hyp=max(abs(y(jj+1)-y(jj)),abs(y(jj)-y(jj-1)));
hym=min(abs(y(jj+1)-y(jj)),abs(y(jj)-y(jj-1)));


hx=5/6*hxp+1/6*hxm+1/9*dx;
hy=5/6*hyp+1/6*hym+1/9*dy;

%hx=1.85*hxp;
%hy=1.85*hyp;

%hx=1.0*hxp;
%hy=1.0*hyp;

%dxlocal=abs(x(ii+1)-x(ii)); dylocal=abs(y(jj+1)-y(jj));
%crit=0.40;
%if (abs(xt-x(ii))/dxlocal > crit) && (abs(yt-y(jj))/dylocal > crit)
%hx=1.10*hxp;
%hy=1.10*hyp;
%end


%hx=1.5*hxp;
%hy=1.5*hyp;

%hx=max(abs(x(ii+1)-x(ii)),abs(x(ii)-x(ii-2)));
%hy=max(abs(y(jj+1)-y(jj)),abs(y(jj)-y(jj-2)));

nsup=0;
for j=jj-5:jj+5
  dy=abs(yt-y(j));
    for i=ii-5:ii+5
      dx=abs(xt-x(i));
              if dx <= 3*hx/2 & dy <= 3*hy/2
                nsup=nsup+1;
                ix(nsup)=i;
                jy(nsup)=j;
                area(nsup)=0.25*(x(i+1)-x(i-1))*(y(j+1)-y(j-1));
              end
     end
end
