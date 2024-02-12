function [nnd,connec,nel,geom]=getMesh(NXE,NYE,dhx,dhy,Q)


if Q==8



nnd = 0;
k = 0;
for j = 1:NYE
for i=1:NXE
k = k + 1;

n1 = (j-1)*(3*NXE+2)+2*i - 1;n2 = j*(3*NXE+2)+i - NXE - 1;n3 = j*(3*NXE+2)+2*i-1;
n4 = n3 + 1;n5 = n3 + 2;n6 = n2 + 1;n7 = n1 + 2;n8 = n1 + 1;

geom(n1,:) = [(i-1)*dhx (j-1)*dhy ]; %left down corner

geom(n3,:) = [(i-1)*dhx j*dhy];  %left top corner
geom(n2,:) = [(geom(n1,1)+geom(n3,1))/2 (geom(n1,2)+geom(n3,2))/2]; %down middle node
geom(n5,:) = [i*dhx j*dhy  ]; %right top corner
geom(n4,:) = [(geom(n3,1)+ geom(n5,1))/2 (geom(n3,2)+ geom(n5,2))/2]; % right middle node
geom(n7,:) = [i*dhx  (j-1)*dhy  ];  %right down corner
geom(n6,:) = [(geom(n5,1)+ geom(n7,1))/2 (geom(n5,2)+ geom(n7,2))/2];
geom(n8,:) = [(geom(n1,1)+ geom(n7,1))/2 (geom(n1,2)+ geom(n7,2))/2];

nel = k;
nnd = n5;
connec(k,:) = [n1 n8 n7 n6 n5 n4 n3 n2]; %counterclockwise from left down corner
end
end




elseif Q==4

lx=NXE*dhx;
ly=NYE*dhy;


HelpVec = (linspace(0, lx, NXE+1));
nnd = (NXE+1)*(NYE+1);
nel = NXE * NYE;

% Set nodes.
xx = zeros(nnd, 1);
yy = zeros(nnd, 1);
for i = 1 : NYE+1
    yy( (i-1)*(NXE+1)+1 : i*(NXE+1) ) = (i-1)*(ly/NYE);
    xx( (i-1)*(NXE+1)+1 : i*(NXE+1) ) = HelpVec;
end


% Define elements.
connec = zeros(nel, 4);
for i = 1 : NXE
    for j = 1 : NYE
        CurrElem = (j-1)*NXE + i;
        c = (j-1)*(NXE+1)+1 + i-1;
        connec(CurrElem, :) =[c c+1 c+NXE+2 c+NXE+1];
    end
end

geom=[xx yy];

end %if
end %function

