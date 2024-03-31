function z=nu_grid(n,l,fac,unif)
% unif 1 monotone strech
if unif ==1
   n1=n-1;
   eta=0.5*(0:n-1)/n1;
   fun=0.5*(1+tanh(fac*(eta-0.5))/tanh(fac*0.5));
   z=2*l*fun;
elseif unif ==0 % uniform grid
   z=linspace(0,l,n);
elseif unif==2
    z=zeros(n,1);
    if mod(n,2) ==1 % odd no of nodes
       nn=(n-1)/2+1;
       n1=nn-1;
       eta=0.5*(0:nn-1)/n1;
       fun=0.5*(1+tanh(fac*(eta-0.5))/tanh(fac*0.5));
       z(nn:n)=l*fun+l/2;
       for i=1:nn-1
           delta=z(nn+i)-z(nn+i-1);
           z(nn-i)=z(nn-i+1)-delta;
       end
    else % even no of nodes
       disp('not implemented...stopping')
       stop
    end
end
