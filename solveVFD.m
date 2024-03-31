function [U]=solveVFD(PX,PY,EX,EY,f,Vb,hx,hy,alpha);
[m,n]=size(f);
%%% SET BOUNDARY VALUES (DIRICHLET)
nyv=n+2; nxv=m+2;
hy2=hy^2; hx2=hx^2;
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G=zeros(m,n); W=zeros(m,n);
%W=E-alpha*ones(m,n); W=1./W;
G=PX'*f*PY;
for j=1:n
    for i=1:m
        W(i,j)=G(i,j)/(EX(i)+EY(j)-alpha);
    end
end
%W=W.*G;
U=PX*W*PY';
clear G; clear W;
