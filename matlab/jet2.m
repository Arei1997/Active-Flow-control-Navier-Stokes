function J4 = jet2(m)
%JET2   Variant of HSV.
%   JET2(M), a variant of HSV(M), is an M-by-3 matrix containing
%   the colormap used by CONTOUR, SURF and PCOLOR.
%   The colors begin with dark blue, range through shades of
%   blue, cyan, white, yellow and red, and end with dark red.
%   JET, with no arguments, is the same length as the current colormap.
%   Use COLORMAP(JET2).
%
%   See also JET, HSV, HOT, PINK, FLAG, COLORMAP, RGBPLOT.

%   Emanuele Terrile


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:20
J4(i,1)=0;
end
for i=21:26
J4(i,1)=0.1429+J4(i-1,1);
end
for i=27:40
J4(i,1)=1.0;
end
for i=41:53
J4(i,1)=J4(i-1,1)-0.0428;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:6
J4(i,2)=0;
end
for i=7:26
J4(i,2)=0.0476+J4(i-1,2);
end
for i=27
J4(i,2)=1.0;
end
for i=28:47
J4(i,2)=J4(i-1,2)-0.0476;
end
for i=48:53
J4(i,2)=0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=14:27
J4(i,3)=1;
end
for i=13:-1:1
J4(i,3)=J4(i+1,3)-0.0428;
end
for i=28:33
J4(i,3)=-0.1429+J4(i-1,3);
end
for i=34:53
J4(i,3)=0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

