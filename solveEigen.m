function kfactors=solveEigen(K,KG,prescribedDof,w,L,E,E1,E2,v12,v21,h,v,nLambda,gdof,activeDof,connec,geom,nnd,isort,nodvalue,heavynodes,tipnodes,tc,globxyGp,sigmanGP)
K(prescribedDof,:)=[];
K(:,prescribedDof)=[];
KG(prescribedDof,:)=[];
KG(:,prescribedDof)=[];

% perform eigenproblem
[Mshape,L1] = eigs(K,KG,nLambda,0);
L1 = diag(L1);

if isort==1
    kfactors=12*(1-v^2)*L1*((2*w)^2)/pi^2/E/h^3;
elseif isort==2
    kfactors=12*(1-v12*v21)*L1*((2*w)^2)/pi^2/E1/h^3;
    %kfactors=12*(1-v12*v21)*L1*((2*w)^2)/pi^2/E2/h^3;
%     kfactors=12*(1-v12*v21)*L1/pi^2/E2/h^3;
end
% kfactors=real(kfactors);
GModShape=zeros(gdof,nLambda);

n=length(activeDof);
for i=1:n
    GModShape(activeDof(i),:)=Mshape(i,:);
end



nHG=numel(heavynodes);

for k=1:nLambda
    L_shape=GModShape(:,k);
    Scale=min([0.15*L,0.15*w])*(1/max(abs(L_shape(:,1))));
    Wnodes=zeros(nnd,1);
    
    for i=1:nnd
        switch nodvalue(i)
            case 0
                Wnodes(i)=L_shape(3*i-2);
            case 1
                F=find(heavynodes==i);
                Wnodes(i)=L_shape(3*i-2)+L_shape(3*nnd+3*F-2);
            case 2
                F=find(tipnodes==i);
                ddF=3*nnd+3*nHG+9*F;
                Wnodes(i)=L_shape(3*i-2)+L_shape(ddF-8);
        end
    end
    
%     Wnodes=real(Wnodes);
%     [B,TF] = rmoutliers(Wnodes,'mean');
%     Wnodes(TF==1)=0;
    geom3=real([geom 2*Scale*Wnodes]);
    figure
    p=patch('Faces', connec, 'Vertices', geom3,'Marker','o','MarkerSize',2,'FaceVertexCData',10*Wnodes);
    p.FaceColor = 'interp';
    %     shading interp
    view(3)
    axis([0 2*w 0 2*L -1 1])
    axis equal
    title([' Mode Shape number = ',num2str(k),'----- Kfactor= ',num2str(round(kfactors(k),4))]);
end


% plot stress contour at gaus points

plotstress(globxyGp,sigmanGP);

end %function



