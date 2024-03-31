function [epsilon,b,ix,jy,area,nsup,hx,hy]=precompute_eps_norkpm(xx,yy,xlag,ylag,dS)
deltax=xx(2)-xx(1); deltay=yy(2)-yy(1);
nls=length(xlag);
ni=length(xx);
nj=length(yy);
ix=zeros(16,nls);
jy=zeros(16,nls);
area=zeros(16,nls);
nsup=zeros(nls,1);
hx=zeros(nls,1);
hy=zeros(nls,1);
for ilag=1:nls
       xll=xlag(ilag);
       yll=ylag(ilag);
       ivalue=floor(xll/deltax);
       jvalue=floor(yll/deltay);
       [ii,jj,ar,nsp,h1,h2]=support1(xll,yll,xx,yy,ivalue-4,ivalue+4,jvalue-4,jvalue+4);
%       [ii,jj,ar,nsp,h1,h2]=support(xll,yll,xx,yy,ni,nj);
       for is=1:nsp
           ix(is,ilag)=ii(is);
           jy(is,ilag)=jj(is);
           area(is,ilag)=ar(is);
       end
       nsup(ilag)=nsp;
       hx(ilag)=h1;
       hy(ilag)=h2;
end

b=zeros(6,nls);b(1,:)=1;
deltax=xx(2)-xx(1); deltay=yy(2)-yy(1);
epsilon=ones(nls,1)*sqrt(deltax^2+deltay^2)/sqrt(2);
