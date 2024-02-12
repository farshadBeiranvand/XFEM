function [xgp,ygp,wgp,Kinp,B,Bforgauss,elementDof,nDof,globxyGp]=StiffnessIN(Q,p,tipnodes,heavynodes,normEl,heavyEl,...
    heavytipEl,neighbortipEl,tipEl,geom,connec,nnd,nodvalue,x1c,y1c,x2c,y2c,tc,dhx,dhy,h,Ainp,kapa,Bforgauss,globxyGp)

%p=Number of ELEMENT
nHG=numel(heavynodes);
indice=connec(p,:);
unsd=[2*indice-1;2*indice];
elementDof=reshape(unsd,1,[]);


nHel=0;
for i=1:Q
    if nodvalue(indice(i))==1
        nHel=nHel+1;
        F=find(heavynodes==indice(i));
        elementDof=[elementDof,2*nnd+2*F-1,2*nnd+2*F];
    end
end

nTel=0;
for i=1:Q
    if nodvalue(indice(i))==2
        nTel=nTel+1;
        F=find(tipnodes==indice(i));
        ddF=2*nnd+2*nHG+8*F;
        elementDof=[elementDof,ddF-7,ddF-6,ddF-5,ddF-4,ddF-3,ddF-2,ddF-1,ddF];
    end
end

nDof=2*Q+nHel*2+nTel*8;

if numel(tipEl)~=0
    if isempty(intersect(indice,connec(tipEl(1),:)))==0
        xc=x1c;yc=y1c;
    else
        xc=x2c;yc=y2c;
    end
else
    xc=x1c;yc=y1c;
end


[Hn,Fn,nodesenrichvalue]=HFin(Q,geom,indice,nodvalue,x1c,y1c,tc,xc,yc); % H,F for nodes

[xgp,ygp,wgp,xshear,yshear,wshear]=getQuadrature(Q,p,normEl,heavyEl,heavytipEl,neighbortipEl,tipEl,geom,indice,x1c,y1c,tc,xc,yc);
Kinp=zeros(nDof);

for j=1:max(size(wgp))
    xi=xgp(j);
    eta=ygp(j);
    gw=wgp(j);
    
    [shapeFunction]=getShapeFunction(xi,eta,Q);
    [naturalDerivatives,secondNaturalderivatives]=getDerivatives(xi,eta,Q);
    [jac,~,XYDerivatives,secXYDerivatives]=Jacobian(geom(indice,:),naturalDerivatives,secondNaturalderivatives);
    [Xxi,Yeta]=globecoordGp(geom,indice,dhx,dhy,xi,eta);
    
    [Hi,Fi,dxFi,dyFi]=HF_int(Q,nodesenrichvalue,geom,indice,xi,eta,dhx,dhy,x1c,y1c,tc,xc,yc);
    
    [K_int,B]=BmatrixINplane(Q,jac,nodesenrichvalue,XYDerivatives,shapeFunction,h,Ainp,kapa,gw,Hn,Fn,Hi,Fi,dxFi,dyFi);
    
    
    Kinp=Kinp+K_int;
    
    Bforgauss{p,j}=B;
    globxyGp=[globxyGp;Xxi Yeta];
end %gausspoints


end












