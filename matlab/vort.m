function [omega,psi,xo,yo]=vort(u,v,hx,hy)
[mm,nn]=size(u);
ncx=mm-1; ncy=nn-2;
omega=zeros(ncx+1,ncy+1);
for j=1:ncy+1
    for i=1:ncx+1
        omega(i,j)=(v(i+1,j)-v(i,j))/hx-(u(i,j+1)-u(i,j))/hy;
    end
end
for i=1:ncx+1
  xo(i)=(i-1)*hx;
end
  
for j=1:ncy+1
  yo(j)=(j-1)*hy;
end
%%%%%%%%%%%%%%%%
m=ncx-1; n=ncy-1;
A=-(2+2*hy^2/hx^2)*eye(m,m);
A=A+hy^2*diag(diag(eye(m-1,m-1)),1)/hx^2+hy^2*diag(diag(eye(m-1,m-1)),-1)/hx^2;
[P,D]=eig(A); 
L=diag(D); clear D;
P1=inv(P);
f=zeros(ncx-1,ncy-1);
for j=2:ncy
    for i=2:ncx
        f(i-1,j-1)=-omega(i,j)*hy^2;
    end
end
FT=zeros(ncx-1,ncy-1);
FT=P1*f;
[nn,mm]=size(FT);
UT=zeros(nn,mm);
for i=1:nn
    diago=L(i)*ones(mm,1);
    UT(i,:)=trispecial(diago,FT(i,:),0);
end
psi=P*UT;
clear UT;

