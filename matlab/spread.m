function [array]=spread(wklag,epsilon,nls,ni,nj,xlag,ylag,dS,x,y,b,...
ixg,jyg,areag,nsupg,hxg,hyg)

array=zeros(ni,nj);
   for ilag=1:nls
       xll=xlag(ilag);
       yll=ylag(ilag);
%       [ix,jy,area,nsup,hx,hy]=support(xll,yll,x,y,ni,nj);
         nsup=nsupg(ilag); hx=hxg(ilag); hy=hyg(ilag);
         ix=ixg(1:nsup,ilag); jy=jyg(1:nsup,ilag); area=areag(1:nsup,ilag); 
       deltav=epsilon(ilag)*dS(ilag);
       bb=b(:,ilag);
       for is=1:nsup
          ry=y(jy(is))-yll;
          rx=x(ix(is))-xll;
          rd=delta2d(rx,ry,hx,hy,bb);
%           array(jy(is),ix(is))=array(jy(is),ix(is))+wklag(ilag)*rd*deltav;
           array(ix(is),jy(is))=array(ix(is),jy(is))+wklag(ilag)*rd*deltav;
          end
  end

