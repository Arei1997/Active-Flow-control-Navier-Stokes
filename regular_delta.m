function value=regular_delta(r,h)
ar=abs(r);
if ar > 3/2*h
   value=0;
elseif ar > h/2
   value=(5-3*ar/h-sqrt(-3*(1-(ar/h))^2+1))/(6*h);
else
   value=(1+sqrt(1-3*(ar/h)^2))/(3*h);
end


