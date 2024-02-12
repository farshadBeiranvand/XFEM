function [Hn,Gn,Rn,nodesenrichvalue]=HGR(Q,geom,indice,nodvalue,x1c,y1c,tc,xc,yc)

nodesenrichvalue=zeros(Q,1);% 0 = normal node
Hn=zeros(Q,1);
Rn=zeros(Q,4); % Rn=[R1;R2;R3;R4] R1=[R1(1);R1(2);...:R1(Q)];
Gn=zeros(Q,1);

for z=1:Q
    xx=geom(indice(z),1);
    yy=geom(indice(z),2);
    nodesenrichvalue(z,1)=nodvalue(indice(z));
    [Heavyvalue]=Heavysidefunction(xx,yy,x1c,y1c,tc);
    Hn(z,1)=Heavyvalue;
    
    if nodesenrichvalue(z,1)==2
        [R1f,R2f,R3f,R4f,G1f]=TIPfunction(xx,yy,xc,yc,tc,x1c,y1c);
        Gn(z,1)=G1f;
        Rn(z,1)=R1f;
        Rn(z,2)=R2f;
        Rn(z,3)=R3f;
        Rn(z,4)=R4f;
        
    end % tip functions
end %nodes


end



