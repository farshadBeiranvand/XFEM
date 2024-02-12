function [subnumber,xisubnodes,etasubnodes]=subelementsTip(xc,yc,xld,yld,xlt,ylt,xrd,yrd,xrt,yrt,xcrack,ycrack,xend,yend)

dhx=xrd-xld;
dhy=ylt-yld;
A=dhx*dhy;

%%%%%%%%%%%%%%%%%%%%%%%%%
xrr=[xld xlt xrd xrt xcrack xend];yrr=[yld ylt yrd yrt ycrack yend]; %sort nodes counterclockwise
cx = mean(xrr);cy = mean(yrr);a = atan2(yrr - cy, xrr - cx);
[~, order] = sort(a);xg = xrr(order);yg = yrr(order);
p=[xc yc;xg(1) yg(1);xg(2) yg(2);xg(3) yg(3);xg(4) yg(4);xg(5) yg(5);xg(6) yg(6)];




%%%%%%%%%%%%%%%%%
tris=[1 1 2 3;1 1 3 4;1 1 4 5;1 1 5 6;1 1 6 7;1 1 7 2];

x=zeros(6,4);
y=zeros(6,4);
for i=1:6
x(i,:)=p(tris(i,:),1);
y(i,:)=p(tris(i,:),2);
end

ddd=size(x);n=ddd(1,1);
h=1;
subnumber=0;

for i=1:n
xpoly=x(i,:);ypoly=y(i,:);
a=polyarea(xpoly,ypoly);              %area of triangle
if a/A > 0.0001
subnumber=subnumber+1;
xisubnodes(h,:)=xpoly;
etasubnodes(h,:)=ypoly;
h=h+1;
end
end


%transform(xx,yy) to xi and eta 
n=size(xisubnodes);
i=n(1,1);
j=n(1,2);

xmean=(xld+xrd)/2;
ymean=(yld+ylt)/2;

q=1;
for p=1:i*j
xisubnodes(p)=(xisubnodes(p)-xmean)*2/dhx;
etasubnodes(p)=(etasubnodes(p)-ymean)*2/dhy;
q=q+1;
end


end




