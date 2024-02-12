function [Kbend,Kgeo,Kshear,elementDof,sigma,sigmanGP]=StiffnessOUT(Q,p,nel,nnd,normEl,heavyEl,heavytipEl,neighbortipEl,tipEl,geom,...
    connec,nodvalue,x1c,y1c,x2c,y2c,tc,dhx,dhy,h,kapa,miu,Ds,Ainp,As,Bforgauss,displacements,heavynodes,tipnodes,elementDof_inplane,NG2,NG3,KWfound,GPfound,sigma,sigmanGP)
%p=Number of ELEMENT

nHG=numel(heavynodes);
indice=connec(p,:);
unsd=[3*indice-2;3*indice-1;3*indice];
elementDof=reshape(unsd,1,[]);


nHel=0;
for i=1:Q
if nodvalue(indice(i))==1
nHel=nHel+1;
F=find(heavynodes==indice(i));
elementDof=[elementDof,3*nnd+3*F-2,3*nnd+3*F-1,3*nnd+3*F];
end
end

nTel=0;
for i=1:Q
if nodvalue(indice(i))==2
nTel=nTel+1;
F=find(tipnodes==indice(i));
TT=3*nnd+3*nHG+9*F;
elementDof=[elementDof,TT-8,TT-7,TT-6,TT-5,TT-4,TT-3,TT-2,TT-1,TT];
end
end

nDof=3*Q+nHel*3+nTel*9;

if numel(tipEl)~=0
if isempty(intersect(indice,connec(tipEl(1),:)))==0
xc=x1c;yc=y1c;
else
xc=x2c;yc=y2c;
end
else
    xc=x1c;yc=y1c;
end

[Hn,Gn,Rn,nodesenrichvalue]=HGR(Q,geom,indice,nodvalue,x1c,y1c,tc,xc,yc);

[xgp,ygp,wgp,xshear,yshear,wshear]=getQuadrature(Q,p,normEl,heavyEl,heavytipEl,neighbortipEl,tipEl,geom,indice,x1c,y1c,tc,xc,yc);
Kbend=zeros(nDof);
Kgeo=zeros(nDof);
Kshear=zeros(nDof);

%%%%%%%%%%%%%
%K bend & geo
%%%%%%%%%%%%%
for j=1:max(size(wgp))
xi=xgp(j);
eta=ygp(j);
gw=wgp(j);

[shapeFunction]=getShapeFunction(xi,eta,Q);
[naturalDerivatives,secondNaturalderivatives]=getDerivatives(xi,eta,Q);
[jac,~,XYDerivatives,secXYDerivatives]=Jacobian(geom(indice,:),naturalDerivatives,secondNaturalderivatives);

[Hi,Ri,Gi,dxRi,dyRi,dxGi,dyGi,d2xxGi,d2yyGi,d2xyGi]=...
HGR_int(Q,nodesenrichvalue,geom,indice,xi,eta,dhx,dhy,xc,yc,tc,x1c,y1c);

[NG1]=NETforce(p,j,elementDof_inplane,Bforgauss,displacements,Ainp);
sigma{p,j}=[NG1(1) NG1(2) NG1(4)]; %Sxx Sxy Syy
sigmanGP=[sigmanGP;NG1(1) NG1(2) NG1(4)];
[Kb_int,Kg_int]=Bmatrix(Q,jac,nodesenrichvalue,XYDerivatives,secXYDerivatives,shapeFunction,...
h,miu,Ds,gw,NG1,NG2,NG3,Hn,Rn,Gn,Hi,Ri,Gi,dxRi,dyRi,dxGi,dyGi,d2xxGi,d2yyGi,d2xyGi,KWfound,GPfound);

Kbend=Kbend+Kb_int;
Kgeo=Kgeo+Kg_int;

end


%K shear

for j=1:max(size(wshear))
xi=xshear(j);
eta=yshear(j);
gw=wshear(j);

[shapeFunction]=getShapeFunction(xi,eta,Q);
[naturalDerivatives,secondNaturalderivatives]=getDerivatives(xi,eta,Q);
[jac,~,XYDerivatives,secXYDerivatives]=Jacobian(geom(indice,:),naturalDerivatives,secondNaturalderivatives);

[Hi,Ri,Gi,dxRi,dyRi,dxGi,dyGi,d2xxGi,d2yyGi,d2xyGi]=...
HGR_int(Q,nodesenrichvalue,geom,indice,xi,eta,dhx,dhy,xc,yc,tc,x1c,y1c);

[Ks_int]=BmatrixShear(Q,jac,nodesenrichvalue,XYDerivatives,shapeFunction,h,As,kapa,gw,Hn,Rn,Gn,Hi,Ri,Gi,dxGi,dyGi);

Kshear=Kshear+Ks_int;
end



end



