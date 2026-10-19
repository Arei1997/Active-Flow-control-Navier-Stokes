function [U]=solveUFD(PX,PY,EX,EY,F,Ub,Ubs,hx,hy,alpha);
[m,n]=size(F);
%%% SET BOUNDARY VALUES (DIRICHLET)
hy2=hy^2; hx2=hx^2;
for j=1:n
   F(1,j)=F(1,j)-Ub(4,j+1)/hx2;
   F(m,j)=F(m,j)-Ubs(j)/hx2;
end
G=zeros(m,n); W=zeros(m,n);
%W=E-alpha*ones(m,n); W=1./W;
G=PX'*F*PY;
for j=1:n
    for i=1:m
        W(i,j)=G(i,j)/(EX(i)+EY(j)-alpha);
    end
end
%W=W.*G;
U=PX*W*PY';
clear G; clear W;
