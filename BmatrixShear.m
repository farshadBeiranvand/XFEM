function [Ks_int]=BmatrixShear(Q,jac,nodesenrichvalue,XYDerivatives,shapeFunction,h,As,kapa,gw,Hn,Rn,Gn,Hi,Ri,Gi,dxGi,dyGi)

B1Snorm=zeros(2,3*Q);
for i=1:Q
    
    %%%%%norm
    B1Snorm(1,(i-1)*3+1)=XYDerivatives(i,1);
    B1Snorm(2,(i-1)*3+1)=XYDerivatives(i,2);
    B1Snorm(1,(i-1)*3+2)=shapeFunction(i);
    B1Snorm(2,(i-1)*3+3)=shapeFunction(i);
end
B1S=B1Snorm;


for i=1:Q
    if nodesenrichvalue(i)==1
        suy=(Hi-Hn(i));
        B1SH=suy*B1S(:,3*i-2:3*i);
        B1S=[B1S,B1SH];
    end
end



for k=1:Q
    if nodesenrichvalue(k)==2
        
        R1N=shapeFunction(k)*(Ri(k,1)-Rn(k,1));
        R2N=shapeFunction(k)*(Ri(k,2)-Rn(k,2));
        R3N=shapeFunction(k)*(Ri(k,3)-Rn(k,3));
        R4N=shapeFunction(k)*(Ri(k,4)-Rn(k,4));
        
        G1Nx=XYDerivatives(k,1)*(Gi(k)-Gn(k))+(dxGi(k))*shapeFunction(k);
        G1Ny=XYDerivatives(k,2)*(Gi(k)-Gn(k))+(dyGi(k))*shapeFunction(k);
        
        
        B1St=[...
            G1Nx R1N 0 R2N 0 R3N 0 R4N 0
            G1Ny 0 R1N 0 R2N 0 R3N 0 R4N]; %G2=G3=G4=0 --> wt2 wt3 wt4 is removed
        
        
        B1S=[B1S,B1St];
    end
    
    
    
    
    dtj=det(jac);
    
    Ks_int=kapa*B1S'*As*B1S*gw*dtj;
    
end


