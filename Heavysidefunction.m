function [Heavyvalue]=Heavysidefunction(xx,yy,x1c,y1c,tc)
ln=(xx-x1c)/cosd(tc);
if yy>y1c+ln*sind(tc)
Heavyvalue=1;
elseif yy<y1c+ln*sind(tc)
Heavyvalue=-1;
else
Heavyvalue=0;
end %if
end


