function [b]=assemble_rkpm(x,y,xlag,ylag,dS,ixg,jyg,areag,nsupg,hxg,hyg)
% find the correction polynomial
nls=length(xlag);
ni=length(x);nj=length(y);
msize=6;
   [powx,powy]=powmat;
   b=zeros(msize,nls);
   for ilag=1:nls
        xll=xlag(ilag);   yll=ylag(ilag);
%        xll=0.5*(xlag(ilag)+xlag(ilag+1));
%        yll=0.5*(ylag(ilag)+ylag(ilag+1));
%       [ix,jy,area,nsup,hx,hy]=support(xll,yll,x,y,ni,nj);
       nsup=nsupg(ilag);
       ix=ixg(1:nsup,ilag); jy=jyg(1:nsup,ilag); area=areag(1:nsup,ilag); 
       hx=hxg(ilag); hy=hyg(ilag);
       M=zeros(msize,msize);
       [M,H]=moments(x,y,xll,yll,ix,jy,area,nsup,hx,hy,powx,powy);
       dm=det(M); %disp(nsup);
       if dm==0
          b(1:msize,ilag)=0;
          b(1,ilag)=1;
       else
          p=zeros(msize,1); p(1)=1;
          LM= chol(M,'lower');
          yM=LM\p;
          b(:,ilag)=LM'\yM;
%          b(:,ilag)=M\p;
          b(:,ilag)=H*b(:,ilag);
       end
    end
end
