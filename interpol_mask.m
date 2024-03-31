function [valout]=interpol_mask(array,nls,xlag,ylag,dS,x,y,b,...
ixg,jyg,areag,nsupg,hxg,hyg,rc,xc,yc,iwhat)
hx=x(2)-x(1); hy=y(2)-y(1);
valout=zeros(nls,1);
if iwhat==1
   [n,m]=size(array);
   tmp=zeros(n,m);
   for i=2:n-1
       xp=x(i-1);
       for j=2:m-1
           yp=y(j-1);
           dist=sqrt((xp-xc)^2+(yp-yc)^2);
           if dist < rc
              tmp(i,j)=array(i,j);
           end
       end   
   end 
   lphi=laplaphi(tmp,hx,hy);
   for ilag=1:nls
       xll=xlag(ilag);
       yll=ylag(ilag);
        nsup=nsupg(ilag);
        ix=ixg(1:nsup,ilag); jy=jyg(1:nsup,ilag); area=areag(1:nsup,ilag); 
        hx=hxg(ilag); hy=hyg(ilag); 
       bb=b(:,ilag);
         for is=1:nsup
	  xp=x(ix(is));yp=y(jy(is));
	  ry=yp-yll;
          rx=xp-xll;
          rd=delta2d(rx,ry,hx,hy,bb);
          dist=sqrt((xp-xc)^2+(yp-yc)^2);
          if dist > rc
             rd=0;
          end
          valout(ilag)=valout(ilag)+lphi(ix(is)+1,jy(is)+1)*rd*area(is);
          end
    end
else
   for ilag=1:nls
       xll=xlag(ilag);
       yll=ylag(ilag);
        nsup=nsupg(ilag);
        ix=ixg(1:nsup,ilag); jy=jyg(1:nsup,ilag); area=areag(1:nsup,ilag); 
        hx=hxg(ilag); hy=hyg(ilag); 
       bb=b(:,ilag);
         for is=1:nsup
	  xp=x(ix(is));yp=y(jy(is));
	  ry=yp-yll;
          rx=xp-xll;
          rd=delta2d(rx,ry,hx,hy,bb);
          dist=sqrt((xp-xc)^2+(yp-yc)^2);
          if dist > rc
             rd=0;
          end
          valout(ilag)=valout(ilag)+array(ix(is),jy(is))*rd*area(is);
          end
    end
end
