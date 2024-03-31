function f=sol_P_stag(P,V,L,EE,f,hx,hy)
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
hy2=hy^2; hx2=hx^2;
%
FT=P'*f;
UT=zeros(n,m);
ff=zeros(m,1); a=zeros(m,1); alpha=zeros(m,1); beta=zeros(m,1);
for i=1:n-1
    a=ones(m,1)*L(i);
    a(1)=a(1)+1/hy^2;
    a(m)=a(m)+1/hy^2;
%    keyboard
    [alpha,beta]=specialLU(a,1/hy^2);
    ff=FT(i,:);
    uu=specialbid(alpha,beta,ff);
    UT(i,:)=uu;
end
i=n;
ff=FT(i,:);
u1=V'*ff'; 
for j=1:m; 
    u1(j)=u1(j)/EE(j);
end
uu=V*u1;
UT(i,:)=uu;
f=P*UT;
