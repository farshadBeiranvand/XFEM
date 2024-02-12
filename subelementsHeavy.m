function [subnumber,xisubnodes,etasubnodes]=subelementsHeavy(xld,yld,xlt,ylt,xrd,yrd,xrt,yrt,xcrack,ycrack,xend,yend)

p=[xld yld;xlt ylt;xrd yrd;xrt yrt;xcrack ycrack;xend yend];
dhx=xrd-xld;
dhy=ylt-yld;
A=dhx*dhy;

if ycrack==yld || ycrack==ylt   
first=1; %horizontal
else
first=0; %vertical
end
if yend==yld || yend==ylt
second=1;
else
second=0;
end


h=1;
subnumber=0;
if (first==1 && second==0) || (first==0 && second==1) %triangular subelements
del=delaunay(p);
m=size(del);n=m(1,1);
for i=1:n
x(i,:)=p(del(i,:),1);
y(i,:)=p(del(i,:),2);
end
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




%#########################################
else %quadrilateral subelements

x=[xld xlt xrd xrt xcrack xend];y=[yld ylt yrd yrt ycrack yend]; %sort nodes counterclockwise
cx = mean(x);cy = mean(y);a = atan2(y - cy, x - cx);
[~, order] = sort(a);x = x(order);y = y(order);
for i=1:6
if x(i)==xcrack && y(i)==ycrack
s=i;
t=6-i;
end
end

z=1;
x2=zeros(1,6);
y2=zeros(1,6);
for k=1:t+1
x2(k)=x(s);
y2(k)=y(s);
s=s+1;
end
for j=t+2:6
x2(j)=x(z);
y2(j)=y(z);
z=z+1;
end
x3=[x2(1) x2(2) x2(3) x2(4);x2(4) x2(5) x2(6) x2(1)];
y3=[y2(1) y2(2) y2(3) y2(4);y2(4) y2(5) y2(6) y2(1)];
for i=1:2
xpoly=x3(i,:);ypoly=y3(i,:);
a=polyarea(xpoly,ypoly);              %area of triangle
if a/A > 0.0001
subnumber=subnumber+1;
xisubnodes(h,:)=xpoly;
etasubnodes(h,:)=ypoly;
h=h+1;
end
end

end
%transform(xx,yy) to xi and eta 
nn=size(xisubnodes);
i=nn(1,1);j=nn(1,2);

xmean=(xld+xrd)/2;
ymean=(yld+ylt)/2;

q=1;
for ps=1:i*j
xisubnodes(ps)=(xisubnodes(ps)-xmean)*2/dhx;
etasubnodes(ps)=(etasubnodes(ps)-ymean)*2/dhy;
q=q+1;
end






end




