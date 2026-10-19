      function out1=precondfun(vec)
       global nls ni nj lx ly xc yc rc
       global xlag ylag dS flag x y b stop
       out=zeros(nls,1);
       vec=vec(:);
  for kl=1:nls
  if abs(flag(kl))>stop	  
   window_centered=zeros(nj,ni);
   xkl=xlag(kl);   ykl=ylag(kl);
   [ix,jy,Asup,nsupkl,hx,hy]=support(xkl,ykl);
   bb=b(:,kl);
   sum_windows=0;
     for is=1:nsupkl
      ry=y(jy(is))-ykl;  rx=x(ix(is))-xkl;
      window_centered(jy(is),ix(is))=delta2d(rx,ry,hx,hy,bb);
      dsum=Asup(is) * ( window_centered(jy(is),ix(is)) )^2;
      sum_windows=sum_windows+dsum;
     end
  diag_precond(kl)=(dS(kl)*sum_windows);
  else
  diag_precond(kl)=1;
  end
  end
diag_precond=diag_precond(:);
%out1=spdiags(ones(nls,1)*diag_precond.^(1/2),0,nls,nls);

%out=ones(1,nls)./diag_precond;
out1=vec./(diag_precond.^(1));
out1=out1(:);

