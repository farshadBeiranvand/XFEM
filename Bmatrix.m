function [Kb_int,Kg_int]=Bmatrix(Q,jac,nodesenrichvalue,XYDerivatives,secXYDerivatives,shapeFunction,...
    h,miu,Ds,gw,NG1,NG2,NG3,Hn,Rn,Gn,Hi,Ri,Gi,dxRi,dyRi,dxGi,dyGi,d2xxGi,d2yyGi,d2xyGi,KWfound,GPfound)

B1Gnorm=zeros(2,3*Q);
B2Snorm=zeros(3,3*Q);
B2Gnorm=zeros(3,3*Q);

B1Fnorm=zeros(2,3*Q); %foundation matrixes
B2Fnorm=zeros(2,3*Q);
for i=1:Q %cycle for nodes to define B matrices containing shapeFunction and naturalDerivatives
    
    %%%normal
    B2Snorm(1,(i-1)*3+2)=XYDerivatives(i,1);
    B2Snorm(3,(i-1)*3+3)=XYDerivatives(i,1);
    B2Snorm(2,(i-1)*3+3)=XYDerivatives(i,2);
    B2Snorm(3,(i-1)*3+2)=XYDerivatives(i,2);
    
    B1Gnorm(1,(i-1)*3+1)=XYDerivatives(i,1);
    B1Gnorm(2,(i-1)*3+1)=XYDerivatives(i,2);
    
    B2Gnorm(1,(i-1)*3+1)=secXYDerivatives(i,1);
    B2Gnorm(2,(i-1)*3+1)=secXYDerivatives(i,3);
    B2Gnorm(3,(i-1)*3+1)=2*secXYDerivatives(i,2);
    
    
    B1Fnorm(1,(i-1)*3+1)=shapeFunction(i);
    B2Fnorm(1,(i-1)*3+1)=XYDerivatives(i,1);
    B2Fnorm(2,(i-1)*3+1)=XYDerivatives(i,2);
end

B2S=B2Snorm;
B1G=B1Gnorm;
B2G=B2Gnorm;

B1F=B1Fnorm;
B2F=B2Fnorm;

%%%%% BH
for i=1:Q
    if nodesenrichvalue(i)==1
        suy=(Hi-Hn(i));
        
        B2SH=suy*B2S(:,3*i-2:3*i);
        B1GH=suy*B1G(:,3*i-2:3*i);
        B2GH=suy*B2G(:,3*i-2:3*i);
        B1FH=zeros(2,3);
        B2FH=zeros(2,3);
        
        B2S=[B2S,B2SH];
        B1G=[B1G,B1GH];
        B2G=[B2G,B2GH];
        B1F=[B1F,B1FH];%foundation
        B2F=[B2F,B2FH];
        
        
        
    end
end



for k=1:Q
    if nodesenrichvalue(k)==2
        
        R1Nx=XYDerivatives(k,1)*(Ri(k,1)-Rn(k,1))+(dxRi(k,1))*shapeFunction(k);
        R1Ny=XYDerivatives(k,2)*(Ri(k,1)-Rn(k,1))+(dyRi(k,1))*shapeFunction(k);
        
        R2Nx=XYDerivatives(k,1)*(Ri(k,2)-Rn(k,2))+(dxRi(k,2))*shapeFunction(k);
        R2Ny=XYDerivatives(k,2)*(Ri(k,2)-Rn(k,2))+(dyRi(k,2))*shapeFunction(k);
        
        R3Nx=XYDerivatives(k,1)*(Ri(k,3)-Rn(k,3))+(dxRi(k,3))*shapeFunction(k);
        R3Ny=XYDerivatives(k,2)*(Ri(k,3)-Rn(k,3))+(dyRi(k,3))*shapeFunction(k);
        
        R4Nx=XYDerivatives(k,1)*(Ri(k,4)-Rn(k,4))+(dxRi(k,4))*shapeFunction(k);
        R4Ny=XYDerivatives(k,2)*(Ri(k,4)-Rn(k,4))+(dyRi(k,4))*shapeFunction(k);
        
        G1Nx=XYDerivatives(k,1)*(Gi(k)-Gn(k))+(dxGi(k))*shapeFunction(k);
        G1Ny=XYDerivatives(k,2)*(Gi(k)-Gn(k))+(dyGi(k))*shapeFunction(k);
        
        G1Nxx=secXYDerivatives(k,1)*(Gi(k)-Gn(k))+2*((dxGi(k))*XYDerivatives(k,1))+(d2xxGi(k))*shapeFunction(k);
        G1Nyy=secXYDerivatives(k,3)*(Gi(k)-Gn(k))+2*((dyGi(k))*XYDerivatives(k,2))+(d2yyGi(k))*shapeFunction(k);
        G1Nxy=2*(secXYDerivatives(k,2)*(Gi(k)-Gn(k))+XYDerivatives(k,1)*(dyGi(k))+XYDerivatives(i,2)*(dxGi(k))+...
            (d2xyGi(k))*shapeFunction(k));
        
        B2St=[...
            0 R1Nx 0    R2Nx 0    R3Nx 0    R4Nx 0
            0 0    R1Ny 0    R2Ny 0    R3Ny 0    R4Ny
            0 R1Ny R1Nx R2Nx R2Ny R3Nx R3Ny R4Nx R4Ny];
        
        
        B1Gt=[...
            G1Nx 0 0 0 0 0 0 0 0
            G1Ny 0 0 0 0 0 0 0 0]; % G2=G3=G4=0 --> wt2 wt3 wt4 is removed
        
        B2Gt=[...
            G1Nxx 0 0 0 0 0 0 0 0 
            G1Nyy 0 0 0 0 0 0 0 0 
            G1Nxy 0 0 0 0 0 0 0 0 ]; % G2=G3=G4=0 --> wt2 wt3 wt4 is removed
        
        B1Ft=zeros(2,9);%foundation
        B2Ft=zeros(2,9);
        
        
        B2S=[B2S,B2St];
        B1G=[B1G,B1Gt];
        B2G=[B2G,B2Gt];
        B1F=[B1F,B1Ft];
        B2F=[B2F,B2Ft];
        
        
    end %if
end %for

%Kbending & KG for element

dtj=det(jac);
Kfound=-KWfound*(B1F'*B1F)*gw*dtj-GPfound*(B2F'*B2F)*gw*dtj;
Kb_int=B2S'*Ds*B2S*gw*dtj+Kfound;

Kg_int=B1G'*NG1*B1G*gw*dtj+miu*B2G'*NG2*B2G*gw*dtj+2*miu*B2G'*NG3*B2G*gw*dtj;

end



