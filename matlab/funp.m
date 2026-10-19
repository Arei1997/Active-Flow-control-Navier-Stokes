function [u,f]=funp(x,y,alpha)
u=y^2*(y/0.3e1-0.1e1/0.2e1)*x^2*(x/0.3e1-0.1e1/0.2e1)+cos(pi*x)*cos(pi*y);
f=0.2e1/0.3e1*x^3*y-x^2*y-x^3/0.3e1+x^2/0.2e1-0.2e1*cos(pi*x)*cos(pi*y)* ...
pi^2+0.2e1/0.3e1*y^3*x-y^3/0.3e1-y^2*x+y^2/0.2e1-alpha*u;

