
function [xxIntRef, yyIntRef, wwIntRef] = IntPoints2DCrackTipElem(xxElem,yyElem, xcrack,ycrack,xend,yend, xc, yc, nQ1D);
%xc,yc,xld,yld,xlt,ylt,xrd,yrd,xrt,yrt,xcrack,ycrack,xend,yend

ffElem=zeros(1,4);
for i=1:4
a=xxElem(i);b=yyElem(i);
s=[xend-xcrack a-xcrack;yend-ycrack b-ycrack];
ffElem(i)=sign(det(s));
end




SignVect = sign(ffElem);
Count = 0;
if SignVect(1) ~= SignVect(2)
    Count = Count + 1; 
    CutSegm(Count) = 1;
end 
if SignVect(2) ~= SignVect(3)
    Count = Count + 1; 
    CutSegm(Count) = 2;
end 
if SignVect(3) ~= SignVect(4)
    Count = Count + 1; 
    CutSegm(Count) = 3;
end 
if SignVect(4) ~= SignVect(1)
    Count = Count + 1; 
    CutSegm(Count) = 4;
end 


% Construct the six integration triangles ...
% ... number the nodes.
Count = 1;
for i = 1 : 4
    HelpNodes(i+Count-1) = i;
    if Count <= 2
        if CutSegm(Count) == i
            HelpNodes(i+Count) = 4+Count;
            Count = Count + 1;
        end
    end
end


[xcRef, ycRef, CaseIntPointInElem] =ProjectRealToRefElem(xxElem, yyElem, xc, yc);
[xcrackRef,ycrackRef, CaseIntPointInElem] =ProjectRealToRefElem(xxElem,yyElem,xcrack,ycrack);
[xendRef,yendRef, CaseIntPointInElem] =ProjectRealToRefElem(xxElem,yyElem,xend,yend);

if xendRef>xcrackRef
xxTri = [[-1  1 1 -1],xendRef, xcrackRef, xcRef];
yyTri = [[-1 -1 1  1],yendRef, ycrackRef, ycRef];
else
xxTri = [[-1  1 1 -1],xcrackRef, xendRef, xcRef];
yyTri = [[-1 -1 1  1],ycrackRef,yendRef, ycRef];
end



for i = 1 : 6
    MeshTri(i, 1) = 7;
    MeshTri(i, 2) = HelpNodes(i);
    if i == 6
        MeshTri(i, 3) = HelpNodes(1);
    else
        MeshTri(i, 3) = HelpNodes(i+1);
    end
end


% Set integration points in each triangle according to almost polar integration.
xxIntRef = zeros(1, 6*nQ1D*nQ1D);
yyIntRef = zeros(1, 6*nQ1D*nQ1D);
wwIntRef = zeros(1, 6*nQ1D*nQ1D);
[xxIntRefTri, yyIntRefTri, wwIntRefTri] = IntPoints2DRefElemTri(nQ1D);
for i = 1 : 6
    NodesTri = MeshTri(i, :);
    xxElemTri = xxTri(NodesTri);
    yyElemTri = yyTri(NodesTri);
    [xxInt, yyInt, wwInt] = IntPoints2DRealElemTri(xxElemTri, yyElemTri, ...
        xxIntRefTri, yyIntRefTri, wwIntRefTri, nQ1D*nQ1D);
    nQ2D = nQ1D * nQ1D;
    xxIntRef((i-1)*nQ2D+1 : i*nQ2D) = xxInt;
    yyIntRef((i-1)*nQ2D+1 : i*nQ2D) = yyInt;
    wwIntRef((i-1)*nQ2D+1 : i*nQ2D) = wwInt;
end

end

	