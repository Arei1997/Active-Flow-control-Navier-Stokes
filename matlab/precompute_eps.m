function [epsilon,b,ix,jy,area,nsup,hx,hy]=precompute_eps(xx,yy,xlag,ylag,dS)

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

%bound_cond ='D';
%deltn=2*pi*rc/nls;%min([min(diff(x)),min(diff(y))])^2;%

%b=zeros(6,nls);b(1,:)=1;

b=assemble_rkpm(xx,yy,xlag,ylag,dS,ix,jy,area,nsup,hx,hy);

%fun=uact/deltat;
[flag]=ones(nls,1);

bic_b=ones(nls,1);
%------------------------------- matrix formulation
M=buildmat(nls,ni,nj,xlag,ylag,dS,flag,xx,yy,b,ix,jy,area,nsup,hx,hy);
eta=M\bic_b;
epsilon=eta(:);
save M.mat M;
