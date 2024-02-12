% Copyright (c) 2010, Thomas-Peter Fries, RWTH Aachen University
function [xxInt, yyInt, wwInt] = IntPoints2DRealElemQuad(xxElem, yyElem, ...
    xxIntRef, yyIntRef, wwIntRef, nQxy)

% Get nQxy integration points and weights in a real quad-element with  
% element nodes at (xxElem, yyElem).
%
% Example call:
% nQxy = 11;
% [xxIntRef, yyIntRef, wwIntRef] = IntPoints2DRefElemQuad(nQxy);
% [xxInt, yyInt, wwInt] = IntPoints2DRealElemQuad([4 -5 -3 2], ...
%     [3 4 -4 -3], xxIntRef, yyIntRef, wwIntRef, nQxy*nQxy)

[N, dNdx, dNdy, xxInt, yyInt, wwInt] = ShapeFctsStrdFEM(...
    xxElem, yyElem, xxIntRef, yyIntRef, wwIntRef, nQxy);

% % Plot situation.
% reset(cla), reset(clf), hold on
% patch(xxElem, yyElem, 'y')
% a = plot(xxInt, yyInt, 'ko');
% set(a, 'MarkerSize', 8, 'MarkerFaceColor', 'k')
% axis equal
