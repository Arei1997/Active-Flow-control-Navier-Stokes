function [alpha,beta]=specialLU(a,extra)
%%%%%
% given 
%           a(1)   extra 0.................0
%           extra   a(2) extra 0...........0
%  A=       0      extra  a(3) extra 0.....0
%
%           0......................extra a(n)
% returns alpha and beta such that A=B*C
% with
%           1       0............................0
%           beta(2) 1 0..........................0
%  B=       0      beta(3) 1 0 0.................0
%
%           0..........................0 beta(n) 1
%
%           alpha(1) extra 0.....................0
%           0       alpha(2) extra 0 ............0
%  C=       0           0    alpha(3) extra 0....0
%
%           0..........................  0 alpha(n) 
% NOTE that beta(1) is used to store "extra"
%
n=length(a); 
alpha=zeros(n,1);
beta=zeros(n,1);
alpha(1)=a(1);
for k=2:n
    beta(k)=extra/alpha(k-1);
    alpha(k)=a(k)-extra^2/alpha(k-1);
end
beta(1)=extra;
