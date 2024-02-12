function[xcrack,ycrack,xend,yend,xxElem,yyElem]=firstendpoint(Q,geom,indice,x1c,y1c,tc)

if Q==8
xld=geom(indice(1),1);yld=geom(indice(1),2);%left down corner
xrd=geom(indice(3),1);yrd=geom(indice(3),2);
xrt=geom(indice(5),1);yrt=geom(indice(5),2);% right top corner
xlt=geom(indice(7),1);ylt=geom(indice(7),2);
elseif Q==4
xld=geom(indice(1),1);yld=geom(indice(1),2);
xrd=geom(indice(2),1);yrd=geom(indice(2),2);
xrt=geom(indice(3),1);yrt=geom(indice(3),2);
xlt=geom(indice(4),1);ylt=geom(indice(4),2);
end



xxElem=[xld,xrd,xrt,xlt];
yyElem=[yld,yrd,yrt,ylt];




dhx=xrd-xld;
dhy=ylt-yld;
p=3*(dhx^2+dhy^2);
x1c=x1c-p*cosd(tc);
y1c=y1c-p*sind(tc);

ln=(xld-x1c)/cosd(tc);
ylocal=y1c+ln*sind(tc);

if tc==90
xcrack=x1c;ycrack=yld;xend=x1c;yend=ylt;

else

if ylocal>ylt %1
ycrack=ylt; xcrack=(y1c-ycrack)*cotd(-1*tc)+x1c; %(x,ycrack are coordinates of crackpoint that comeinto the element)
if (xrt-xcrack)*tand(-1*tc)<dhy
xend=xrt;yend=ycrack-(xrt-xcrack)*tand(-1*tc);   %(x,yend are coordinates of crackpoint that leave the element)
elseif (xrt-xcrack)*tand(-1*tc)==dhy
xend=xrd;yend=yrd;
else
yend=yrd;xend=xcrack+dhy*cotd(-1*tc);
end


elseif ylocal<ylt && ylocal>yld     %2         
xcrack=xld; ycrack=ylocal;
if tc>0
if(ylt-ycrack)*cotd(tc)<dhx
yend=ylt;xend=xlt+(ylt-ycrack)*cotd(tc);
elseif (ylt-ycrack)*cotd(tc)>dhx
xend=xrt;yend=ycrack+dhx*tand(tc);
else
yend=yrt;xend=xrt;
end
end

if tc==0
xend=xcrack+dhx;
yend=ycrack;
end


if tc<0
if(ycrack-yld)*cotd(-1*tc)<dhx
yend=yld;xend=xcrack+(ycrack-yld)*cotd(-1*tc);
elseif (ycrack-yld)*cotd(-1*tc)>dhx
xend=xrd;yend=ycrack-dhx*tand(-1*tc);
else
yend=yrd;xend=xrd;
end
end



elseif ylocal<yld %3
ycrack=yld; xcrack=(ycrack-y1c)*cotd(tc)+x1c;
if (xrd-xcrack)*tand(tc)<dhy
xend=xrd;yend=ycrack+(xrd-xcrack)*tand(tc);
elseif (xrd-xcrack)*tand(tc)>dhy
yend=yrt;xend=xcrack+dhy*cotd(tc);
else
xend=xrt;yend=yrt;
end


elseif ylocal==yld %4
ycrack=yld;xcrack=xld;
if (xrd-xcrack)*tand(tc)<dhy
xend=xrd;yend=ycrack+(xrd-xcrack)*tand(tc);
elseif (xrd-xcrack)*tand(tc)>dhy
yend=yrt;xend=xcrack+dhy*cotd(tc);
else
xend=xrt;yend=yrt;
end


elseif ylocal==ylt %5
ycrack=ylt;xcrack=xlt;
if (xrt-xcrack)*tand(-1*tc)<dhy
xend=xrt;yend=ycrack-(xrt-xcrack)*tand(-1*tc);   %(x,yend are coordinates of crackpoint that leave the element)
elseif(xrt-xcrack)*tand(-1*tc)==dhy
xend=xrd;yend=yrd;
else
yend=yrd;xend=xcrack+dhy*cotd(-1*tc);
end

end
end





end


