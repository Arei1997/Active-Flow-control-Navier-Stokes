%-------------------------------------------------------------------------
function [xlag,ylag,theta,dS]=initialize_lag(xx,yy,rc,xc,yc,irkpm,actime)
% lagrangian parameters
dxmin=min(diff(xx));
dymin=min(diff(yy));
diagm=sqrt(dxmin^2+dymin^2);
nls=floor(0.95*2*pi*rc/diagm);
if irkpm==0
    disp('precomputed nls ')
    disp(nls)
    nlsnew=input('if you want to change nls enter the new value, otherwise enter 0 ')
    if nlsnew >= 0 
       nls=nlsnew
    end
end

[theta,xlag,ylag,dS]=setupdisk(nls,xc,yc,rc,actime);

%xlag=xlag+rc*cos(2*pi*omega*actime);
%ylag=ylag+rc*sin(2*pi*omega*actime);

%-------------------------------------------------------------------------
function [theta,xlag,ylag,dS]=setupdisk(nls,xc,yc,rc,actime)
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

