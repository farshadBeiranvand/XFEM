function [F1,F2,F3,F4,dxF1,dxF2,dxF3,dxF4,dyF1,dyF2,dyF3,dyF4]=TIPfunctionINplane(xxi,yyi,xc,yc,tc,x1c,y1c)
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


% thStar must be in the range from -pi to +pi.


% Define the four branch functions and derivatives.
%F1
F1=sqrt(rr)*sin(thStar/2); 
F1_X1=-0.5*rr^-0.5*sin(thStar/2);
F1_X2=0.5*rr^-0.5*cos(thStar/2);
dxF1=F1_X1*cos(Alpha)-F1_X2*sin(Alpha);
dyF1=F1_X1*sin(Alpha)+F1_X2*cos(Alpha);

%F2
F2=sqrt(rr)*cos(thStar/2); 
F2_X1=0.5*rr^-0.5*cos(thStar/2);
F2_X2=0.5*rr^-0.5*sin(thStar/2); 
dxF2=F2_X1*cos(Alpha)-F2_X2*sin(Alpha);
dyF2=F2_X1*sin(Alpha)+F2_X2*cos(Alpha); 

%F3
F3=sqrt(rr)*sin(thStar/2)*sin(thStar);
F3_X1=-0.5*rr^-0.5*sin(1.5*thStar)*sin(thStar);
F3_X2=0.5*rr^-0.5*(sin(thStar/2)+sin(1.5*thStar)*cos(thStar));
dxF3=F3_X1*cos(Alpha)-F3_X2*sin(Alpha);
dyF3=F3_X1*sin(Alpha)+F3_X2*cos(Alpha); 

%F4
F4=sqrt(rr)*cos(thStar/2)*sin(thStar);
F4_X1=-0.5*rr^-0.5*cos(1.5*thStar)*sin(thStar);
F4_X2=0.5*rr^-0.5*(cos(thStar/2)+cos(1.5*thStar)*cos(thStar));
dxF4=F4_X1*cos(Alpha)-F4_X2*sin(Alpha);
dyF4=F4_X1*sin(Alpha)+F4_X2*cos(Alpha);

end
