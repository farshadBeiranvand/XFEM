function [Hi,Ri,Gi,dxRi,dyRi,dxGi,dyGi,d2xxGi,d2yyGi,d2xyGi]=...
HGR_int(Q,nodesenrichvalue,geom,indice,xi,eta,dhx,dhy,xc,yc,tc,x1c,y1c)


Ri=zeros(Q,4);
Gi=zeros(Q,1);

dxRi=zeros(Q,4);
dyRi=zeros(Q,4);
dxGi=zeros(Q,1);
dyGi=zeros(Q,1);

d2xxGi=zeros(Q,1);
d2yyGi=zeros(Q,1);
d2xyGi=zeros(Q,1);

if Q==8
xmean=geom(indice(2),1);
ymean=geom(indice(4),2);

elseif Q==4
xmean=(geom(indice(1),1)+geom(indice(2),1))/2;
ymean=(geom(indice(1),2)+geom(indice(4),2))/2;
end


Xxi=xmean+xi*dhx/2;
Yeta=ymean+eta*dhy/2;
[Heavyvalue]=Heavysidefunction(Xxi,Yeta,x1c,y1c,tc);
Hi=Heavyvalue;

%R & G (xi,eta)=Rx Gx
[r1,r2,r3,r4,g1,dxr1,dyr1,dxr2,dyr2,dxr3,dyr3,dxr4,dyr4,dxg1,dyg1,d2xxg1,d2yyg1,d2xyg1]=TIPfunction(Xxi,Yeta,xc,yc,tc,x1c,y1c);

for z=1:Q
if nodesenrichvalue(z,1)==2

Ri(z,1)=r1;
Ri(z,2)=r2;
Ri(z,3)=r3;
Ri(z,4)=r4;
Gi(z,1)=g1;

dxRi(z,1)=dxr1;
dyRi(z,1)=dyr1;
dxRi(z,2)=dxr2;
dyRi(z,2)=dyr2;
dxRi(z,3)=dxr3;
dyRi(z,3)=dyr3;
dxRi(z,4)=dxr4;
dyRi(z,4)=dyr4;

dxGi(z)=dxg1;
dyGi(z)=dyg1;

d2xxGi(z)=d2xxg1;
d2yyGi(z)=d2yyg1;
d2xyGi(z)=d2xyg1;
end
end


end


