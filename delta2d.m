function [rd]=delta2d(rx,ry,hx,hy,b)
deltay=regular_delta(ry,hy);
deltax=regular_delta(rx,hx);
rxx=rx; ryy=ry;
rd=(b(1)+b(2)*rxx+b(3)*ryy+b(4)*rxx*ryy+b(5)*rxx^2+b(6)*ryy^2)*deltax*deltay;
