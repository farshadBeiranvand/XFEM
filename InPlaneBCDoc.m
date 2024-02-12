function [UnodResistant,VnodResistant]=InPlaneBCDoc(NXE,NYE)


if NXE==floor(NXE/2)*2
% Even Mesh
UnodResistant=[(NXE/2)+1,(NXE+1)*NYE+(NXE/2)+1];
else
% Odd Mesh
UnodResistant=[(NXE+1)/2,((NXE+1)/2)+1,(NXE+1)*NYE+(NXE+1)/2,(NXE+1)*NYE+((NXE+1)/2)+1];
end
if NYE==floor(NYE/2)*2
% Even Mesh
VnodResistant=[1+(NXE+1)*NYE/2,(NXE+1)+(NXE+1)*NYE/2]; 
else
% Odd Mesh
VnodResistant=[(NXE+1)*((NYE+1)/2),(NXE+1)*(1+((NYE+1)/2)),(NXE+1)*((NYE+1)/2)-NXE,(NXE+1)*(1+((NYE+1)/2))-NXE]; 
end


end