function [Xxi,Yeta]=globecoordGp(geom,indice,dhx,dhy,xi,eta)
x0=geom(indice(1),1);
y0=geom(indice(1),2);
Xxi=(x0+dhx/2)+xi*dhx/2;
Yeta=(y0+dhy/2)+eta*dhy/2;
end