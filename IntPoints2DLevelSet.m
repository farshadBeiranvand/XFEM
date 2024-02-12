% Copyright (c) 2010, Thomas-Peter Fries, RWTH Aachen University
function [xgp, ygp, wgp] = IntPoints2DLevelSet(xxElem,yyElem,xcrack,ycrack,xend,yend, nQx,nQy)

[xcrackRef,ycrackRef] =ProjectRealToRefElem(xxElem,yyElem,xcrack,ycrack);
[xendRef,yendRef] =ProjectRealToRefElem(xxElem,yyElem,xend,yend);

xA=xcrackRef;
yA=ycrackRef;
xB=xendRef;
yB=yendRef;
x1 =-1; y1 =-1;x2 = 1; y2 =-1;
x3 = 1; y3 = 1;x4 =-1; y4 = 1;

xi=[-1 1 1 -1];
eta=[-1 -1 1 1];
SignVect=zeros(1,4);
for i=1:4
    a=xi(i);b=eta(i);
    s=[xB-xA a-xA;yB-yA b-yA];
    SignVect(i)=sign(det(s));
end

% Set integration points and weigths in the 2D reference element
[xxIntRef, yyIntRef, wwIntRef] = IntPoints2DRefElemQuad(nQx,nQy);
nQxy = nQx* nQy;

if SignVect(1)*SignVect(2)*SignVect(3)*SignVect(4) < 0
    % Treatment of Type-2-elements (one level function value is on the other side than the other three).
    % (-+++, +-++, ++-+, +++-, +---, -+--, --+-, ---+)
    if SignVect(1)~=SignVect(2) && SignVect(1)~=SignVect(3) && SignVect(1)~=SignVect(4)
        pos=1;
    elseif SignVect(2)~=SignVect(1) && SignVect(2)~=SignVect(3) && SignVect(2)~=SignVect(4)
        pos=2;
    elseif SignVect(3)~=SignVect(1) && SignVect(3)~=SignVect(2) && SignVect(3)~=SignVect(4)
        pos=3;
    elseif SignVect(4)~=SignVect(1) && SignVect(4)~=SignVect(2) && SignVect(4)~=SignVect(3)
        pos=4;
    end
    
    xP=xi(pos);
    yP=eta(pos);
    
    xC = 0.5*(xA+xB);
    yC = 0.5*(yA+yB);
    xD=.5*(xP+xC);
    yD=.5*(yP+yC);
    xF=.5*(xP+xB);
    yF=.5*(yP+yB);
    xE=.5*(xP+xA);
    yE=.5*(yP+yA);
    
    
    
    if pos==1
        [xxInt1, yyInt1, wwInt1] = IntPoints2DRealElemQuad([x4 xA xC x3], [y4 yA yC y3], xxIntRef, yyIntRef, wwIntRef, nQxy);
        [xxInt2, yyInt2, wwInt2] = IntPoints2DRealElemQuad([x3 xC xB x2], [y3 yC yB y2], xxIntRef, yyIntRef, wwIntRef, nQxy);
        [xxInt3, yyInt3, wwInt3] = IntPoints2DRealElemQuad([xA xE xD xC], [yA yE yD yC], xxIntRef, yyIntRef, wwIntRef, nQxy);
        [xxInt4, yyInt4, wwInt4] = IntPoints2DRealElemQuad([xB xC xD xF], [yB yC yD yF], xxIntRef, yyIntRef, wwIntRef, nQxy);
        [xxInt5, yyInt5, wwInt5] = IntPoints2DRealElemQuad([x1 xF xD xE], [y1 yF yD yE], xxIntRef, yyIntRef, wwIntRef, nQxy);
    elseif pos==2
        [xxInt1, yyInt1, wwInt1] = IntPoints2DRealElemQuad([x1 xA xC x4], [y1 yA yC y4], xxIntRef, yyIntRef, wwIntRef, nQxy);
        [xxInt2, yyInt2, wwInt2] = IntPoints2DRealElemQuad([x4 xC xB x3], [y4 yC yB y3], xxIntRef, yyIntRef, wwIntRef, nQxy);
        [xxInt3, yyInt3, wwInt3] = IntPoints2DRealElemQuad([xA xE xD xC], [yA yE yD yC], xxIntRef, yyIntRef, wwIntRef, nQxy);
        [xxInt4, yyInt4, wwInt4] = IntPoints2DRealElemQuad([xB xC xD xF], [yB yC yD yF], xxIntRef, yyIntRef, wwIntRef, nQxy);
        [xxInt5, yyInt5, wwInt5] = IntPoints2DRealElemQuad([x2 xF xD xE], [y2 yF yD yE], xxIntRef, yyIntRef, wwIntRef, nQxy);
    elseif pos==3
        [xxInt1, yyInt1, wwInt1] = IntPoints2DRealElemQuad([x4 xA xC x1], [y4 yA yC y1], xxIntRef, yyIntRef, wwIntRef, nQxy);
        [xxInt2, yyInt2, wwInt2] = IntPoints2DRealElemQuad([x1 xC xB x2], [y1 yC yB y2], xxIntRef, yyIntRef, wwIntRef, nQxy);
        [xxInt3, yyInt3, wwInt3] = IntPoints2DRealElemQuad([xA xE xD xC], [yA yE yD yC], xxIntRef, yyIntRef, wwIntRef, nQxy);
        [xxInt4, yyInt4, wwInt4] = IntPoints2DRealElemQuad([xB xC xD xF], [yB yC yD yF], xxIntRef, yyIntRef, wwIntRef, nQxy);
        [xxInt5, yyInt5, wwInt5] = IntPoints2DRealElemQuad([x3 xF xD xE], [y3 yF yD yE], xxIntRef, yyIntRef, wwIntRef, nQxy);
    elseif pos==4
        [xxInt1, yyInt1, wwInt1] = IntPoints2DRealElemQuad([x1 xA xC x2], [y1 yA yC y2], xxIntRef, yyIntRef, wwIntRef, nQxy);
        [xxInt2, yyInt2, wwInt2] = IntPoints2DRealElemQuad([x2 xC xB x3], [y2 yC yB y3], xxIntRef, yyIntRef, wwIntRef, nQxy);
        [xxInt3, yyInt3, wwInt3] = IntPoints2DRealElemQuad([xA xE xD xC], [yA yE yD yC], xxIntRef, yyIntRef, wwIntRef, nQxy);
        [xxInt4, yyInt4, wwInt4] = IntPoints2DRealElemQuad([xB xC xD xF], [yB yC yD yF], xxIntRef, yyIntRef, wwIntRef, nQxy);
        [xxInt5, yyInt5, wwInt5] = IntPoints2DRealElemQuad([x4 xF xD xE], [y4 yF yD yE], xxIntRef, yyIntRef, wwIntRef, nQxy);
    end
    
    xgp = [xxInt1, xxInt2, xxInt3, xxInt4, xxInt5];
    ygp = [yyInt1, yyInt2, yyInt3, yyInt4, yyInt5];
    wgp = [wwInt1, wwInt2, wwInt3, wwInt4, wwInt5];
    %     plot(xxInt1, yyInt1, 'k*')
    %     plot(xxInt2, yyInt2, 'k*')
    %     plot(xxInt3, yyInt3, 'k*')
    %     plot(xxInt4, yyInt4, 'k*')
    %     plot(xxInt5, yyInt5, 'k*')
    
    
else
    % Treatment of Type-3-elements (two level function values are on the other side than the other two)
    % (+--+, ++--, +-+-, -++-, --++, -+-+)
    
    if SignVect(1)==SignVect(2) && SignVect(3)==SignVect(4)
        [xxInt1, yyInt1, wwInt1] = IntPoints2DRealElemQuad([x1 x2 xB xA], [y1 y2 yB yA], xxIntRef, yyIntRef, wwIntRef, nQxy);
        [xxInt2, yyInt2, wwInt2] = IntPoints2DRealElemQuad([xA xB x3 x4], [yA yB y3 y4], xxIntRef, yyIntRef, wwIntRef, nQxy);
    elseif SignVect(1)==SignVect(4) && SignVect(2)==SignVect(3)
        [xxInt1, yyInt1, wwInt1] = IntPoints2DRealElemQuad([x4 x1 xA xB ], [y4 y1 yA yB ], xxIntRef, yyIntRef, wwIntRef, nQxy);
        [xxInt2, yyInt2, wwInt2] = IntPoints2DRealElemQuad([x3 xB xA x2 ], [y3 yB yA y2 ], xxIntRef, yyIntRef, wwIntRef, nQxy);
    end
    
    
    
    xgp = [xxInt1, xxInt2];
    ygp = [yyInt1, yyInt2];
    wgp = [wwInt1, wwInt2];
    
%         patch([x1 xA xB x4], [y1 yA yB y4], 'b')
    %     patch([xB xA x2 x3], [yB yA y2 y3], 'b')
    %
    %     plot(xxInt1, yyInt1, 'k*')
    %     plot(xxInt2, yyInt2, 'k*')
    
    
end


end