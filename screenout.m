function [dummy]=screenout(uold,u,vold,v,actime,hx,hy,ncx,ncy,cd,cl)
nxu=ncx+1; nyu=ncy+2;
nxv=ncx+2; nyv=ncy+1;
disp('/****** at time ***/');
actime
difu=abs(uold-u); difv=abs(vold-v);
difumax=max(max(difu));
difvmax=max(max(difv));
divmax=divergence_check(u,v,hx,hy,ncx,ncy);
vinlet=sum(u(1,2:nyu-1));
voutlet=sum(u(nxu,2:nyu-1));
disp('max div');
disp(divmax)
disp('max dif in u');
disp(difumax)
disp('max dif in v');
disp(difvmax)
disp('CD');
disp(cd)
disp('CL');
disp(cl)
%disp('sum vel inlet');
%disp(vinlet)
%disp('sum vel outlet');
%disp(voutlet)
dummy=0;

