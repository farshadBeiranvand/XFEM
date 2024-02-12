function [K_int,B]=BmatrixINplane(Q,jac,nodesenrichvalue,XYDerivatives,shapeFunction,h,Ainp,kapa,gw,Hn,Fn,Hi,Fi,dxFi,dyFi)

Bnorm=zeros(3,2*Q);

%%%%% B standard
for i=1:Q
    Bnorm(1,2*i-1)=XYDerivatives(i,1);
    Bnorm(3,2*i)=XYDerivatives(i,1);
    Bnorm(2,2*i)=XYDerivatives(i,2);
    Bnorm(3,2*i-1)=XYDerivatives(i,2);
end
B=Bnorm;


%%%%% BH
for i=1:Q
    if nodesenrichvalue(i)==1
        suy=-abs((Hi-Hn(i)));
        BH=suy*Bnorm(:,2*i-1:2*i);
        B=[B,BH];
    end
end

%%%%% Bt

for k=1:Q
    if nodesenrichvalue(k)==2
        
        F1Nx=XYDerivatives(k,1)*(Fi(k,1)-Fn(k,1))+dxFi(k,1)*shapeFunction(k);
        F1Ny=XYDerivatives(k,2)*(Fi(k,1)-Fn(k,1))+dyFi(k,1)*shapeFunction(k);
        
        F2Nx=XYDerivatives(k,1)*(Fi(k,2)-Fn(k,2))+dxFi(k,2)*shapeFunction(k);
        F2Ny=XYDerivatives(k,2)*(Fi(k,2)-Fn(k,2))+dyFi(k,2)*shapeFunction(k);
        
        F3Nx=XYDerivatives(k,1)*(Fi(k,3)-Fn(k,3))+dxFi(k,3)*shapeFunction(k);
        F3Ny=XYDerivatives(k,2)*(Fi(k,3)-Fn(k,3))+dyFi(k,3)*shapeFunction(k);
        
        F4Nx=XYDerivatives(k,1)*(Fi(k,4)-Fn(k,4))+dxFi(k,4)*shapeFunction(k);
        F4Ny=XYDerivatives(k,2)*(Fi(k,4)-Fn(k,4))+dyFi(k,4)*shapeFunction(k);
        
        
        Bt=[F1Nx,0,F2Nx,0,F3Nx,0,F4Nx,0
            0,F1Ny,0,F2Ny,0,F3Ny,0,F4Ny
            F1Ny,F1Nx,F2Ny,F2Nx,F3Ny,F3Nx,F4Ny,F4Nx ];
        B=[B,Bt];
    end
end



dtj=det(jac);
K_int=B'*Ainp*B*gw*dtj;
end


