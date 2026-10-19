function [theta,xlag,ylag,dS]=setupdisk_fourier(nls,xc,yc,rc)
dtheta=2*pi/nls; 
% s counterclockwise 
for i=1:nls
    theta(i)=dtheta*(i-1)+dtheta; %+(0.5-rand(1))*0.25*dtheta;
    arg=theta(i);
    xlag(i)=xc+cos(arg)*rc;
    ylag(i)=yc+sin(arg)*rc;
end
thp=0.5*(theta(2)+theta(1));
thm=0.5*(theta(1)+theta(nls)-2*pi);
dS(1)=rc*(thp-thm);
for i=2:nls-1
    thp=0.5*(theta(i+1)+theta(i));
    thm=0.5*(theta(i)+theta(i-1));
    dS(i)=rc*(thp-thm);
end
thp=0.5*(theta(1)+2*pi+theta(nls));
thm=0.5*(theta(nls)+theta(nls-1));
dS(nls)=rc*(thp-thm);
%dS(1:nls)=2*pi*rc/nls;
