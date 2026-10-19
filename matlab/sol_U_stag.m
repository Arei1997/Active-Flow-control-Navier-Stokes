function f=sol_U_stag(P,L,f,Ub,hx,hy)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                               %
%               North=3                                                         %
%           ¡---------------------¡                                             %
%   West=4  ¡                     ¡ East=2                                      %
%           ¡                     ¡                                             %
%           -----------------------                                             %
%               South=1                                                         %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
[n,m]=size(f);
nyu=m+2; nxu=n+2;
hy2=hy^2; hx2=hx^2;
%%% SET BOUNDARY VALUES (DIRICHLET)
% west boundary
i=1;
j=1;
f(i,j)=f(i,j)-Ub(4,j+1)/hx2-2*Ub(1,i+1)/hy2;
for j=2:nyu-3
   f(i,j)=f(i,j)-Ub(4,j+1)/hx2;
end
j=nyu-2;
f(i,j)=f(i,j)-Ub(4,j+1)/hx2-2*Ub(3,i+1)/hy2;
%
% east boundary
i=nxu-2;
j=1;
f(i,j)=f(i,j)-Ub(2,j+1)/hx2-2*Ub(1,i+1)/hy2;
for j=2:nyu-3
    f(i,j)=f(i,j)-Ub(2,j+1)/hx2;
end
j=nyu-2;
f(i,j)=f(i,j)-2*Ub(3,i+1)/hy2-Ub(2,j+1)/hx2;
% south boundary
j=1;
for i=2:nxu-3
    f(i,j)=f(i,j)-2*Ub(1,i+1)/hy2;
end
% north boundary
j=nyu-2;
for i=2:nxu-3
    f(i,j)=f(i,j)-2*Ub(3,i+1)/hy2;
end
%
FT=P'*f;
UT=zeros(n,m);
ff=zeros(m,1); a=zeros(m,1); alpha=zeros(m,1); beta=zeros(m,1);
for i=1:n
    a=ones(m,1)*L(i);
    a(1)=a(1)-1/hy^2;
    a(m)=a(m)-1/hy^2;
    [alpha,beta]=specialLU(a,1/hy^2);
    ff=FT(i,:);
    uu=specialbid(alpha,beta,ff);
    UT(i,:)=uu;
end
f=P*UT;
