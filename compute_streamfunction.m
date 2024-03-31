function rhs=compute_rhs_psi(lx,ncx,ly,ncy,u,v)
	hx=lx/ncx; hy=ly/ncy;
	nxu=size(u,1);nxv=size(v,1);
	nyu=size(u,2);nyv=size(v,2);
        %streamfunction	
        for i=2:nxu-1
         for j=2:nyu-2
         dudy(i-1,j-1)=(u(i,j)-u(i,j-1))/hy;
         end
        end
        for i=2:nxv-2
         for j=2:nyv-1
         dvdx(i-1,j-1)=(v(i,j)-v(i-1,j))/hx;
         end
        end
	rhs= -dudy+dvdx;

