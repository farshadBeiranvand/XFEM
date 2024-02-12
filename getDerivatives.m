function[naturalDerivatives,secondNaturalderivatives]=getDerivatives(x,y,Q)

if Q==8
naturalDerivatives=[1/4*(1-y)*(2*x+y) 1/4*(1-x)*(x+2*y);
-x*(-y+1) -1/2*(1+x)*(1-x);
1/4*(1-y)*(2*x-y) 1/4*(1+x)*(2*y-x);
1/2*(1+y)*(1-y) -y*(x+1);
1/4*(1+y)*(2*x+y) 1/4*(1+x)*(x+2*y);
-x*(y+1) 1/2*(1+x)*(1-x);
1/4*(1+y)*(2*x-y) 1/4*(1-x)*(2*y-x);
-1/2*(1+y)*(1-y) -y*(-x+1)];

secondNaturalderivatives=[(1-y)/2 1/4*(-2*y-2*x+1) (1-x)/2; %[xx;xy;yy]
y-1 x 0;
(1-y)/2 -x/2+y/2-1/4 (x+1)/2;
0 -y -x-1;%4
(y+1)/2 x/2+y/2+1/4 (x+1)/2;
-y-1 -x 0;
(y+1)/2 x/2-y/2-1/4 (1-x)/2;
0 y x-1];

elseif Q==4
naturalDerivatives=[-0.25*(1-y) -0.25*(1-x);
+0.25*(1-y) -0.25*(1+x);
+0.25*(1+y) +0.25*(1+x);
-0.25*(1+y) +0.25*(1-x)];

secondNaturalderivatives=[0 0 .025;%[xx;xy;yy]
0 0 -.25;
0 0 .25;
0 0 -.25];

end

end