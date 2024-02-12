function [Hi,Fi,dxFi,dyFi]=HF_int(Q,nodesenrichvalue,geom,indice,xi,eta,dhx,dhy,x1c,y1c,tc,xc,yc)

Fi=zeros(Q,4);
dxFi=zeros(Q,4);
dyFi=zeros(Q,4);

if Q==8
xmean=geom(indice(2),1);
ymean=geom(indice(4),2);

elseif Q==4
xmean=(geom(indice(1),1)+geom(indice(2),1))/2;
ymean=(geom(indice(1),2)+geom(indice(4),2))/2;
end


Xxi=xmean+xi*dhx/2;Yeta=ymean+eta*dhy/2;
[Heavyvalue]=Heavysidefunction(Xxi,Yeta,x1c,y1c,tc);
Hi=Heavyvalue;

%Fi =F_int_point (xi,eta)
[f1,f2,f3,f4,dxf1,dxf2,dxf3,dxf4,dyf1,dyf2,dyf3,dyf4]=TIPfunctionINplane(Xxi,Yeta,xc,yc,tc,x1c,y1c);


for z=1:Q
if nodesenrichvalue(z,1)==2
Fi(z,1)=f1;
Fi(z,2)=f2;
Fi(z,3)=f3;
Fi(z,4)=f4;

dxFi(z,1)=dxf1;
dxFi(z,2)=dxf2;
dxFi(z,3)=dxf3;
dxFi(z,4)=dxf4;

dyFi(z,1)=dyf1;
dyFi(z,2)=dyf2;
dyFi(z,3)=dyf3;
dyFi(z,4)=dyf4;
end
end



end



