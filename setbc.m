function [u,v]=setbc(us,vs,Ub,Vb,Ubs,Vbs)
[nx,ny]=size(us);
nxu=nx+2; nyu=ny+2;
[nx,ny]=size(vs);
nxv=nx+2; nyv=ny+2;
u=zeros(nxu,nyu); v=zeros(nxv,nyv);
u(2:nxu-1,2:nyu-1)=us;
v(2:nxv-1,2:nyv-1)=vs;
for j=1:nyu
    u(1,j)=Ub(4,j);
    u(nxu,j)=Ubs(j);
end
for i=2:nxu
    u(i,1)=u(i,2);
    u(i,nyu)=u(i,nyu-1);
end
for i=1:nxv
    v(i,1)=Vb(1,i);
    v(i,nyv)=Vb(3,i);
end
for j=2:nyv-1
    v(1,j)=2*Vb(4,j)-v(2,j);
    v(nxv,j)=2*Vbs(j)-v(nxv-1,j);
%    v(nxv,j)=2*Vb(2,j)-v(nxv-1,j);
end
