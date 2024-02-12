function [Hn,Fn,nodesenrichvalue]=HFin(Q,geom,indice,nodvalue,x1c,y1c,tc,xc,yc)


nodesenrichvalue=zeros(Q,1);% 0 = normal node
Hn=zeros(Q,1);
Fn=zeros(Q,4); % Fn=[F1;F2;F3;F4] F1=[F1(1);F1(2);...:F1(Q)];





for z=1:Q %cycle for node
xx=geom(indice(z),1);
yy=geom(indice(z),2);
nodesenrichvalue(z,1)=nodvalue(indice(z));
%%%%%%%%%%%%%%Heavyside function
[Heavyvalue]=Heavysidefunction(xx,yy,x1c,y1c,tc);
Hn(z,1)=Heavyvalue;

%%%%%%tip functions
if nodesenrichvalue(z,1)==2

[f1n,f2n,f3n,f4n]=TIPfunctionINplane(xx,yy,xc,yc,tc,x1c,y1c);

Fn(z,1)=f1n;
Fn(z,2)=f2n;
Fn(z,3)=f3n;
Fn(z,4)=f4n;

end
end



end




