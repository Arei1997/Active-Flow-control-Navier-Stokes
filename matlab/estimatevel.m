clear all
rc=0.5; xc=3.003; yc=6.93; Re=30;
irkpm=1; 
load xu; load xv;
load yu; load yv;
load xp; load yp;
load uu; load vv;
hx=xp(2)-xp(1); hy=yp(2)-yp(1);
load pp;
[xlag,ylag,theta,dS]=initialize_lag(xu,yu,rc,xc,yc,irkpm,0);
[epsilonx,bcoefx,ix_x,jy_x,area_x,nsup_x,hx_x,hy_x]=precompute_eps(xu(2:length(xu)-1),yu(2:length(yu)-1),xlag,ylag,dS);
[epsilony,bcoefy,ix_y,jy_y,area_y,nsup_y,hx_y,hy_y]=precompute_eps(xv(2:length(xv)-1),yv(2:length(yv)-1),xlag,ylag,dS);
[epsilonp,bcoefp,ix_p,jy_p,area_p,nsup_p,hx_p,hy_p]=precompute_eps(xp(2:length(xp)-1),yp(2:length(yp)-1),xlag,ylag,dS);
tmp=u(2:length(xu)-1,2:length(yu)-1);
ulagr=interpolante(tmp,xu(2:length(xu)-1),yu(2:length(yu)-1),xlag,ylag,...
                         dS,bcoefx,ix_x,jy_x,area_x,nsup_x,hx_x,hy_x);
clear tmp;
tmp=v(2:length(xv)-1,2:length(yv)-1);
vlagr=interpolante(tmp,xv(2:length(xv)-1),yv(2:length(yv)-1),xlag,ylag,...
                          dS,bcoefy,ix_y,jy_y,area_y,nsup_y,hx_y,hy_y);



