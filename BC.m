
function [prescribedDof,activeDof,fixedNodeW]=BC(BCtype,gdof,geom,nodvalue,nnd,heavynodes,tipnodes)

%% 1 Boundary Condition
xx=geom(:,1);
yy=geom(:,2);

bcR=BCtype(1);
bcD=BCtype(2);
bcL=BCtype(3);
bcU=BCtype(4);

RightNodes=find(xx==max(xx));
DownNodes=find(yy==min(yy));
LeftNodes=find(xx==min(xx));
UpNodes=find(yy==max(yy));

fixedNodeW =[];
fixedNodeTX =[];
fixedNodeTY =[];



switch bcR
	case 's'
		fixedNodeW=[fixedNodeW' RightNodes'];
		fixedNodeTY=[fixedNodeTY RightNodes'];
	case 'c'
		fixedNodeW=[fixedNodeW' RightNodes'];
		fixedNodeTY=[fixedNodeTY' RightNodes'];
		fixedNodeTX=[fixedNodeTX' RightNodes'];
	
end

switch bcD
	case 's'
		fixedNodeW=[fixedNodeW DownNodes'];
		fixedNodeTX=[fixedNodeTX DownNodes'];
	case 'c'
		fixedNodeW=[fixedNodeW DownNodes'];
		fixedNodeTY=[fixedNodeTY DownNodes'];
		fixedNodeTX=[fixedNodeTX DownNodes'];

end


switch bcL
	case 's'
		fixedNodeW=[fixedNodeW LeftNodes'];
		fixedNodeTY=[fixedNodeTY LeftNodes'];
	case 'c'
		fixedNodeW=[fixedNodeW LeftNodes'];
		fixedNodeTY=[fixedNodeTY LeftNodes'];
		fixedNodeTX=[fixedNodeTX LeftNodes'];

end



switch bcU
	case 's'
		fixedNodeW=[fixedNodeW UpNodes'];
		fixedNodeTX=[fixedNodeTX UpNodes'];
	case 'c'
		fixedNodeW=[fixedNodeW UpNodes'];
		fixedNodeTY=[fixedNodeTY UpNodes'];
		fixedNodeTX=[fixedNodeTX UpNodes'];
	
end





z=unique(fixedNodeW);
x=unique(fixedNodeTX);
y=unique(fixedNodeTY);


nHG=numel(heavynodes);

%%W
perW=[];
for i=1:numel(z)
    nod=z(i);
    switch nodvalue(nod)
        case 0
            perW=[perW,3*nod-2];
        case 1
            F=find(heavynodes==nod);
            perW=[perW,3*nod-2,3*nnd+3*F-2];
        case 2
            F=find(tipnodes==nod);
            ddF=3*nnd+3*nHG+9*F;
            perW=[perW,3*nod-2,ddF-8];
    end
end

%%TX
perX=[];
for j=1:numel(x)
    nod=x(j);
    switch nodvalue(nod)
        case 0
            perX=[perX,3*nod-1];
        case 1
            F=find(heavynodes==nod);
            perX=[perX,3*nod-1,3*nnd+3*F-1];
        case 2
            F=find(tipnodes==nod);
            ddF=3*nnd+3*nHG+9*F;
            perX=[perX,3*nod-1,ddF-7,ddF-5,ddF-3,ddF-1];
    end
end


%%TY
perY=[];
for k=1:numel(y)
    nod=y(k);
    switch nodvalue(nod)
        case 0
            perY=[perY,3*nod];
        case 1
            F=find(heavynodes==nod);
            perY=[perY,3*nod,3*nnd+3*F];
        case 2
            F=find(tipnodes==nod);
            ddF=3*nnd+3*nHG+9*F;
            perY=[perY,3*nod,ddF-6,ddF-4,ddF-2,ddF];
    end
end




prescribedDof=[perW,perX,perY];
prescribedDof=unique(prescribedDof);
activeDof=setdiff((1:gdof)',(prescribedDof));


end



