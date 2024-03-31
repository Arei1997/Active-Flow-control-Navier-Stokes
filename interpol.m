function [valout]=interpol(array,nls,ni,nj,xlag,ylag,dS,x,y,b,...
ixg,jyg,areag,nsupg,hxg,hyg)
valout=zeros(nls,1);
   for ilag=1:nls
       xll=xlag(ilag);
       yll=ylag(ilag);
%       [ix,jy,area,nsup,hx,hy]=support(xll,yll,x,y,ni,nj);
        nsup=nsupg(ilag);
        ix=ixg(1:nsup,ilag); jy=jyg(1:nsup,ilag); area=areag(1:nsup,ilag); 
        hx=hxg(ilag); hy=hyg(ilag); 
       bb=b(:,ilag);
         for is=1:nsup
	  xp=x(ix(is));yp=y(jy(is));
	  ry=yp-yll;
          rx=xp-xll;
          rd=delta2d(rx,ry,hx,hy,bb);
%          valout(ilag)=valout(ilag)+array(jy(is),ix(is))*rd*area(is);
          valout(ilag)=valout(ilag)+array(ix(is),jy(is))*rd*area(is);
          end
    end
