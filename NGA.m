function [As,Ds,Ainp,C,NG2,NG3]=NGA(kapa,h,isort,E,v,E1,E2,v12,v21,G12,G13,G23,Nxxnl,Nxynl,Nyynl)


if isort==1

G=E/(2*(1+v));
A11=E*h/(1-v^2);A22=A11;A12=v*E*h/(1-v^2);A44=E/(1-v^2)*0.5*(1-v)*h;A55=A44;A66=A44;
D11=E*h^3/(12*(1-v^2));D22=D11;D12=v*E*h^3/(12*(1-v^2));D66=G*h^3/12;
C=[E/(1-v^2) v*E/(1-v^2) 0;v*E/(1-v^2) E/(1-v^2) 0;0 0 E/2/(1+v)];

elseif isort==2
A11=E1*h/(1-v12*v21);A22=E2*h/(1-v12*v21);A12=v12*E2*h/(1-v12*v21);A44=G23*h;A55=G13*h;A66=G12*h;
D11=E1*h^3/(12*(1-v12*v21));D22=E2*h^3/(12*(1-v12*v21));D12=v12*E2*h^3/(12*(1-v12*v21));D66=G12*h^3/12;
C=zeros(3);
end

As=[A44 0;0 A55];
Ainp=[A11 A12 0;A12 A22 0;0 0 A66];
Ds=[D11 D12 0;D12 D22 0;0 0 D66];
NG2=[Nxxnl Nxxnl 0;Nyynl Nyynl 0;Nxynl Nxynl 0];
NG3=[0 -Nxxnl 1/2*Nxynl;-Nyynl 0 1/2*Nxynl;-1/2*Nxynl -1/2*Nxynl (Nxxnl+Nyynl)/4];

end

