function [M,testdiag]=buildmat(nls,ni,nj,xlag,ylag,dS,flag,x,y,b,ixg,jyg,areag,nsupg,hxg,hyg)
    
 M=zeros(nls,nls);
   test=zeros(ni,nj);
 for ll=1:nls
  window_lcentered=zeros(ni,nj);
  Fl=flag(ll);
  xll=xlag(ll);   yll=ylag(ll);
%  [ixl,jyl,Asupl,nsupl,hx,hy]=support(xll,yll,x,y,ni,nj);
  nsupl=nsupg(ll);
  ixl=ixg(1:nsupl,ll); jyl=jyg(1:nsupl,ll); Asupl=areag(1:nsupl,ll); 
  hx=hxg(ll); hy=hyg(ll);
   bb=b(:,ll);
      for is=1:nsupl
        ry=y(jyl(is))-yll;  rx=x(ixl(is))-xll;
        window_lcentered(ixl(is),jyl(is))=delta2d(rx,ry,hx,hy,bb);
      end

%  testdiag(ll)=4*hx*hy/(dS(kk));
   testdiag(ll)=dS(ll)*buildmat_analytic(bb,hx,hy);
   nblobs=5;
   numblobs=-nblobs:1:nblobs;
   for m=1:2*nblobs+1
   kkvec(m)=[per(ll+numblobs(m),nls)];
   end

   for kk=kkvec
%   for kk=1:nls
   window_kcentered=zeros(ni,nj);
    Fk=flag(kk);
    xkk=xlag(kk);   ykk=ylag(kk);
%    [ixk,jyk,Asupk,nsupk,hx,hy]=support(xkk,ykk,x,y,ni,nj);
     nsupk=nsupg(kk);
     ixk=ixg(1:nsupk,kk); jyk=jyg(1:nsupk,kk); Asupk=areag(1:nsupk,kk); 
     hx=hxg(kk); hy=hyg(kk);
     bb=b(:,kk);
        for is=1:nsupk
          ry=y(jyk(is))-ykk;  rx=x(ixk(is))-xkk;
          window_kcentered(ixk(is),jyk(is))=delta2d(rx,ry,hx,hy,bb);
        end
      sum_windows=0;
      for is=1:nsupl
	dsum=Asupl(is)*window_lcentered(ixl(is),jyl(is))*window_kcentered(ixl(is),jyl(is));
	sum_windows=sum_windows+dsum;
	test(ixl(is),jyl(is))=window_lcentered(ixl(is),jyl(is))*window_kcentered(ixl(is),jyl(is));
      end

    M(ll,kk)=Fk *dS(kk)* sum_windows ;

   end
  end

%    deltav=M\ones(nls,1);
%     LA= chol(AA,'lower');
%     yM=LA\val;
%     deltav=LA'\yM;

%    check_terms

function ind=per(in,ntot)
      ind=in;
      if in>ntot 
      ind=in-ntot;
      elseif in<1
      ind=in+ntot;
      end
