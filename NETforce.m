function [NG1]=NETforce(p,j,elementDof_inplane,Bforgauss,displacements,Ainp)

Be=Bforgauss{p,j};
disp=displacements(elementDof_inplane);
strain=Be*disp;
f=Ainp*strain;

Nxxlo=f(1);
Nyylo=f(2);
Nxylo=f(3);
NG1=-[Nxxlo Nxylo;Nxylo Nyylo];
% NG1=-[0 0;0 1];

end

