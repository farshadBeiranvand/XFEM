function [shapeFunction]=getShapeFunction(x,y,Q)

if Q==8
shapeFunction=[1/4*(1-x)*(1-y)*(-x-y-1);
1/2*(1+x)*(1-y)*(1-x);%2
1/4*(1+x)*(1-y)*(x-y-1);
1/2*(1+x)*(1+y)*(1-y);%4
1/4*(1+x)*(1+y)*(x+y-1);
1/2*(1+x)*(1+y)*(1-x);
1/4*(1-x)*(1+y)*(-x+y-1);%7
1/2*(1-x)*(1+y)*(1-y)];

elseif Q==4
shapeFunction=[0.25*(1-x).*(1-y);
0.25*(1+x).*(1-y);
0.25*(1+x).*(1+y);
0.25*(1-x).*(1+y)];

end

end