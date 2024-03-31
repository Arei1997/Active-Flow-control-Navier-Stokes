function yout=matvec(xin)
      global nls ni nj lx ly xc yc rc  
      global xlag ylag dS flag x y b stop

            yout=interpol(spread(flag,xin));

return

          array=zeros(nj,ni);
	  yout=zeros(nls,1);
	 for ilag=1:nls
            xll=xlag(ilag);
            yll=ylag(ilag);
	    [ix,jy,area,nsup,hx,hy]=support(xll,yll);
            deltav=xin(ilag)*dS(ilag);
            bb=b(:,ilag);
              for is=1:nsup
                  ry=y(jy(is))-yll;
                  rx=x(ix(is))-xll;
                  rd=delta2d(rx,ry,hx,hy,bb);
                  array(jy(is),ix(is))=array(jy(is),ix(is))+flag(ilag)*rd*deltav;
              end
         end

   for ilag=1:nls
       xll=xlag(ilag);
       yll=ylag(ilag);
       [ix,jy,area,nsup,hx,hy]=support(xll,yll);
       bb=b(:,ilag);
       for is=1:nsup
          ry=y(jy(is))-yll;
          rx=x(ix(is))-xll;
          rd=delta2d(rx,ry,hx,hy,bb);
          yout(ilag)=yout(ilag)+array(jy(is),ix(is))*rd*area(is);
       end
   end

