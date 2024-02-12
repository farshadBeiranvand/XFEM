function [R1,R2,R3,R4,G1,dxR1,dyR1,dxR2,dyR2,dxR3,dyR3,dxR4,dyR4,dxG1,dyG1,d2xxG1,d2yyG1,d2xyG1]=TIPfunction(xxi,yyi,xc,yc,tc,x1c,y1c)

tc=tc*pi/180;
if [xc,yc]==[x1c,y1c]
    Alpha=tc+pi;
else
    Alpha=tc;
end
xx = xxi - xc;
yy = yyi - yc;
% Get polar coordinates (rr, thStar) out of (xx, yy).
rr=((xx^2)+(yy^2))^0.5;
if yy>=0 && xx>0
    Tetapr2=(atan(abs(yy)/abs(xx)));
    Tetapr=Tetapr2;
elseif yy>0 && xx<0
    Tetapr2=(atan(abs(yy)/abs(xx)));
    Tetapr=pi-Tetapr2;
elseif yy<=0 && xx<0
    Tetapr2=(atan(abs(yy)/abs(xx)));
    Tetapr=pi+Tetapr2;
elseif yy<0 && xx>0
    Tetapr2=(atan(abs(yy)/abs(xx)));
    Tetapr=2*pi-Tetapr2;
elseif yy>0 && xx==0
    Tetapr=pi/2;
elseif yy<0 && xx==0
    Tetapr=pi+pi/2;
end

teta=double(Tetapr-Alpha);

if teta>pi
    TetaT=teta-2*pi;
elseif teta<-pi
    TetaT=teta+2*pi;
else
    TetaT=teta;
end

thStar=double(TetaT);
% Define the four branch functions and derivatives.
%R1
R1=sqrt(rr)*sin(thStar/2);
R1_X1=-0.5*rr^-0.5*sin(thStar/2);
R1_X2=0.5*rr^-0.5*cos(thStar/2);
dxR1=R1_X1*cos(Alpha)-R1_X2*sin(Alpha);
dyR1=R1_X1*sin(Alpha)+R1_X2*cos(Alpha);

%R2
R2=sqrt(rr)*cos(thStar/2);
R2_X1=0.5*rr^-0.5*cos(thStar/2);
R2_X2=0.5*rr^-0.5*sin(thStar/2);
dxR2=R2_X1*cos(Alpha)-R2_X2*sin(Alpha);
dyR2=R2_X1*sin(Alpha)+R2_X2*cos(Alpha);

%R3
R3=sqrt(rr)*sin(thStar/2)*sin(thStar);
R3_X1=-0.5*rr^-0.5*sin(1.5*thStar)*sin(thStar);
R3_X2=0.5*rr^-0.5*(sin(thStar/2)+sin(1.5*thStar)*cos(thStar));
dxR3=R3_X1*cos(Alpha)-R3_X2*sin(Alpha);
dyR3=R3_X1*sin(Alpha)+R3_X2*cos(Alpha);

%R4
R4=sqrt(rr)*cos(thStar/2)*sin(thStar);
R4_X1=-0.5*rr^-0.5*cos(1.5*thStar)*sin(thStar);
R4_X2=0.5*rr^-0.5*(cos(thStar/2)+cos(1.5*thStar)*cos(thStar));
dxR4=R4_X1*cos(Alpha)-R4_X2*sin(Alpha);
dyR4=R4_X1*sin(Alpha)+R4_X2*cos(Alpha);

%G1
G1=R1;
dxG1=dxR1;
dyG1=dyR1;


d2xxG1=0;
d2yyG1=0;
d2xyG1=0;

end