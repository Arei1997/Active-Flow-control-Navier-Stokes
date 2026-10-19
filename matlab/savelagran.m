function [dummy]=savelagran(u,v,p,xu,yu,xv,yv,xp,yp,xlag,ylag,dS,bcoefx,bcoefy,bcoefp,...
                  ix_x,ix_y,ix_p,jy_x,jy_y,jy_p,area_x,area_y,area_p,...
                  nsup_x,nsup_y,nsup_p,hx_x,hx_y,hx_p,hy_x,hy_y,hy_p,theta)
ulag=interpolante(u,xu,yu,xlag,ylag,dS,bcoefx,ix_x,jy_x,area_x,nsup_x,hx_x,hy_x);
vlag=interpolante(v,xv,yv,xlag,ylag,dS,bcoefy,ix_y,jy_y,area_y,nsup_y,hx_y,hy_y);
plag=interpolante(p,xp,yp,xlag,ylag,dS,bcoefp,ix_p,jy_p,area_p,nsup_p,hx_p,hy_p);
velnorm=ulag.*cos(theta)'+vlag.*sin(theta)';
disp('max normal velocity on Gamma')
disp(max(abs(velnorm)));
disp('max divergence on Gamma')
disp(max(abs(plag)));
save ulag.mat ulag;
save vlag.mat vlag;
save plag.mat plag;
save velnorm.mat velnorm;
dummy=0;
