function [M,H]=moments(x,y,xll,yll,ix,jy,area,nsup,hx,hy,powx,powy)
   msize=6; lend=6;
   H=diag([1,1/hx,1/hy,1/hx/hy,1/hx^2,1/hy^2]);
M=zeros(msize,msize);
for jm=1:lend
    for im=1:lend
        ipow=powx(im,jm);
        jpow=powy(im,jm);
       % loop on i,j on the support mid point rule
        for is=1:nsup
           ii=ix(is);
           ax=x(ii)-xll;
           rd1=regular_delta(ax,hx);
           jj=jy(is);
           ay=y(jj)-yll; 
           rd2=regular_delta(ay,hy);
           aax=ax/hx; aay=ay/hy;
           fun=aax^ipow*aay^jpow*rd1*rd2;
           M(im,jm)=M(im,jm)+area(is)*fun;

%           ax=0.5*(x(ii)+x(ii+1))-xll;
%           rd1=regular_delta(ax,hx);
%           ay=0.5*(y(jj)+y(jj+1))-yll; 
%           rd2=regular_delta(ay,hy);
%           aax=ax/hx; aay=ay/hy;
%           fun=aax^ipow*aay^jpow*rd1*rd2;
%           M(im,jm)=M(im,jm)+0.25*area(is)*fun;
	   
        end
    end
end
