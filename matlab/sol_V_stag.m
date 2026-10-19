function f=sol_V_stag(P,L,f,Vb,hx,hy)
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
nyv=m+2; nxv=n+2;
hy2=hy^2; hx2=hx^2;
%%% SET BOUNDARY VALUES (DIRICHLET)
% west boundary
i=1;
j=1;
f(i,j)=f(i,j)-Vb(1,i+1)/hy2-2*Vb(4,j+1)/hx2;
for j=2:nyv-3
   f(i,j)=f(i,j)-2*Vb(4,j+1)/hx2;
end
j=nyv-2;
f(i,j)=f(i,j)-Vb(3,i+1)/hy2-2*Vb(4,j+1)/hx2;
%
% east boundary
i=nxv-2;
j=1;
f(i,j)=f(i,j)-Vb(1,i+1)/hy2-2*Vb(2,j+1)/hx2;
for j=2:nyv-3
    f(i,j)=f(i,j)-2*Vb(2,j+1)/hx2;
end
j=nyv-2;
f(i,j)=f(i,j)-Vb(3,i+1)/hy2-2*Vb(2,j+1)/hx2;
% south boundary
j=1;
for i=2:nxv-3
    f(i,j)=f(i,j)-Vb(1,i+1)/hy2;
end
% north boundary
j=nyv-2;
for i=2:nxv-3
    f(i,j)=f(i,j)-Vb(3,i+1)/hy2;
end
%
FT=P'*f;
UT=zeros(n,m);
ff=zeros(m,1); a=zeros(m,1); alpha=zeros(m,1); beta=zeros(m,1);
for i=1:n
    a=ones(m,1)*L(i);
%    a(1)=a(1)-1/hx^2;
%    a(m)=a(m)-1/hx^2;
    [alpha,beta]=specialLU(a,1/hy^2);
    ff=FT(i,:);
    uu=specialbid(alpha,beta,ff);
    UT(i,:)=uu;
end
f=P*UT;
