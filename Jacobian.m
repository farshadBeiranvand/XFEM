function [JacobianMatrix,invJacobian,XYDerivatives,secXYDerivatives]=Jacobian(nodeCoordinates,naturalDerivatives,secondNaturalderivatives)

JacobianMatrix=nodeCoordinates'*naturalDerivatives;
invJacobian=inv(JacobianMatrix);
XYDerivatives=naturalDerivatives*invJacobian;


%Second Derivatives 1
Kmatrix=nodeCoordinates'*secondNaturalderivatives;

Ainv=[invJacobian(1,1)^2 2*invJacobian(1,2)*invJacobian(1,1) invJacobian(1,2)^2;
invJacobian(2,1)*invJacobian(1,1) invJacobian(1,1)*invJacobian(2,2)+invJacobian(1,2)*invJacobian(2,1) invJacobian(2,2)*invJacobian(1,2);
invJacobian(2,1)^2 2*invJacobian(2,2)*invJacobian(2,1) invJacobian(2,2)^2];

secXYDerivatives0=(Ainv*secondNaturalderivatives'-Ainv*Kmatrix'*invJacobian*naturalDerivatives')';

%%%%%%%%%%%%%%% 2
j11=invJacobian(1,1);j12=invJacobian(1,2);j21=invJacobian(2,1);j22=invJacobian(2,2);
%k=kesay e=eta
dkk=secondNaturalderivatives(:,1);
dke=secondNaturalderivatives(:,2);
dee=secondNaturalderivatives(:,3);

dxx=j11^2*dkk+2*j11*j12*dke+j12^2*dee;
dxy=j11*j21*dkk+(j11*j22+j12*j21)*dke+j12*j22*dee;
dyy=j21^2*dkk+2*j22*j21*dke+j22^2*dee;
secXYDerivatives=[dxx dxy dyy];

end 