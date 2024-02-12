function plotstress(globxyGp,sigmanGP)
XLocationGp=globxyGp(:,1);
YLocationGp=globxyGp(:,2);

Sxx=sigmanGP(:,1);
Syy=sigmanGP(:,3);
Sxy=sigmanGP(:,2);
% 
% DeltaCxx=max(Sxx)-min(Sxx);
% DeltaCyy=max(Syy)-min(Syy);
% DeltaCxy=max(Sxy)-min(Sxy);

Fxx = scatteredInterpolant(double(XLocationGp),double(YLocationGp),double(Sxx));
Fyy = scatteredInterpolant(double(XLocationGp),double(YLocationGp),double(Syy));
Fxy = scatteredInterpolant(double(XLocationGp),double(YLocationGp),double(Sxy));

[XLocationGpGrid,YLocationGpGrid] = meshgrid(linspace(min(double(XLocationGp')),max(double(XLocationGp'))),linspace(min(double(YLocationGp')),max(double(YLocationGp'))));
%sxx
Sxxgrid = Fxx(XLocationGpGrid,YLocationGpGrid);
figure('Color','white')
contourf(XLocationGpGrid,YLocationGpGrid,Sxxgrid,'LineColor','none','LineStyle','none');
colorbar
hold off
axis equal
title('Sxx Inplane')

%Syy
Syygrid = Fyy(XLocationGpGrid,YLocationGpGrid);
figure('Color','white')
contourf(XLocationGpGrid,YLocationGpGrid,Syygrid,'LineColor','none','LineStyle','none');
colorbar
hold off
axis equal
title('Syy Inplane Without Considing LargeDeformation')

%sxy
Sxygrid = Fxy(XLocationGpGrid,YLocationGpGrid);
figure('Color','white')
contourf(XLocationGpGrid,YLocationGpGrid,Sxygrid,'LineColor','none','LineStyle','none');
colorbar
hold off
axis equal
title('Sxy Inplane Without Considing LargeDeformation')