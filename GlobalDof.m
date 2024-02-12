function [gdof]=GlobalDof(nnd,nodvalue,inout) %1:in 2:out

nH=numel(find(nodvalue==1));
nT=numel(find(nodvalue==2));
if inout==1
gdof=nnd*2+nH*2+nT*8;
elseif inout==2
gdof=nnd*3+nH*3+nT*9;% G2=G3=G4=0 --> wt2,wt3,wt4 is removed
end



end

