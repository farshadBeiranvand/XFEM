function [Nxxnl,Nyynl,Nxynl]=nonlocal_forces(Lfx,Lfy,Nxxlo,Nyylo)
switch Lfx
case 0
Nxxnl=Nxxlo;
case 0.5
Nxxnl=.75*Nxxlo;
case 1
Nxxnl=.5*Nxxlo;
case 1.5
Nxxnl=.3125*Nxxlo;
case 2
Nxxnl=0;
end

switch Lfy
case 0
Nyynl=Nyylo;
case 0.5
Nyynl=.75*Nyylo;
case 1
Nyynl=.5*Nyylo;
case 1.5
Nyynl=.3125*Nyylo;
case 2
Nyynl=0;
end
Nxynl=0;