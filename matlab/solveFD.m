function [U]=solveFD(PX,PY,E,F);
[m,n]=size(F);
G=zeros(m,n); W=zeros(m,n);
G=PX'*F*PY;
G(m,n)=0;
W=E.*G;
U=PX*W*PY';
clear G; clear W;

