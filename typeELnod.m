function [x1c,y1c,x2c,y2c,tipnodes,heavynodes,tipEl,heavyEl,heavytipEl,neighbortipEl,normEl,nodvalue]=typeELnod(geom,nnd,nel,x1c,y1c,lc,tc,neltip,connec,Q,crackedplate)
heavyEl=[];
heavycount=0;

xx=geom(:,1);
yy=geom(:,2);

x2c=x1c+lc*cosd(tc);
y2c=y1c+lc*sind(tc);

c1=[x1c y1c];
c2=[x2c y2c];

n=(c2-c1)/norm(c2-c1)*[0 -1;1 0]; %normal to crack
t1=(c1-c2)/norm(c2-c1); %tangent to c1
t2=(c2-c1)/norm(c2-c1); %tangent to c2

%% Level Set
fn=zeros(nnd,1);%normal level set
ft1=zeros(nnd,1);%ct1 level set
ft2=zeros(nnd,1);%ct2 level set

t1x=t1(1);
t1y=t1(2);
nx=n(1);
ny=n(2);
cx=c1(1);
cy=c1(2);
for i=1:nnd
    x=xx(i);
    y=yy(i);
    if (t1y == 0)
        fn(i)=(cy-y)/ny;
        ft1(i)=(x-cx+fn(i)*nx)/t1x;
    else
        fn(i)=(cx-x-(cy-y)*t1x/t1y)/(nx-ny*t1x/t1y);
        ft1(i)=(y-cy+fn(i)*ny)/t1y;
    end
end

%obtain ft2
t2x=t2(1);
t2y=t2(2);
nx=n(1);
ny=n(2);
cx=c2(1);
cy=c2(2);
for i=1:nnd
    x=xx(i);
    y=yy(i);
    if (t2y == 0)
        fnn=(cy-y)/ny;
        ft2(i)=(x-cx+fnn*nx)/t2x;
    else
        fnn=(cx-x-(cy-y)*t2x/t2y)/(nx-ny*t2x/t2y);
        ft2(i)=(y-cy+fnn*ny)/t2y;
    end
end


Nelem=[];
T1elem=[];
T2elem=[];
Ncount=0;
T1count=0;
T2count=0;
for i=1:nel
    s=connec(i,:);
    if min(sign(fn(s))) ~= max(sign(fn(s)))
        Ncount=Ncount+1;
        Nelem(Ncount)=i;
    end
    
    if min(sign(ft1(s))) ~= max(sign(ft1(s)))
        T1count=T1count+1;
        T1elem(T1count)=i;
    end
    
    if min(sign(ft2(s))) ~= max(sign(ft2(s)))
        T2count=T2count+1;
        T2elem(T2count)=i;
    end
end

CT1=intersect(Nelem,T1elem);
if numel(CT1)>1
    for j=1:numel(CT1)
        nods=connec(CT1(j),:);
        d=sum(fn(nods));
        dis1(j)=d;
    end
index = find(abs(dis1)==min(abs(dis1)));
CT1=CT1(index);
end
    

        
        
CT2=intersect(Nelem,T2elem);
if numel(CT2)>1
    for j=1:numel(CT2)
        nods=connec(CT2(j),:);
        d=sum(fn(nods));
        dis2(j)=d;
    end
index = find(abs(dis2)==min(abs(dis2)));
CT2=CT2(index);
end



if neltip==2
    tipEl=[CT1;CT2];
elseif neltip==1
    tipEl=CT2;
    heavyEl=CT1;
    heavycount=1;
end


tipnodes=unique(connec(tipEl,:));
nlevelnodes=unique(connec(Nelem,:));
ss=size(nlevelnodes);

k=0;
ac=[];
for i=1:ss(1)*ss(2)
    if ft1(nlevelnodes(i))<0 && ft2(nlevelnodes(i))<0
        k=k+1;
        ac(k)=nlevelnodes(i);	 %along crack elements nodes
    end
end

heavynodes=unique(setdiff(ac,tipnodes));
normalnodes=setdiff(setdiff(1:nnd,tipnodes),heavynodes);

%% ElementTYPE
heavytipEl=[];
heavytipcount=0;
neighbortipEl=[];
neighcount=0;
normEl=[];
normcount=0;


for i=1:nel
    nodes=connec(i,:);
    nT=numel(intersect(nodes,tipnodes));
    nH=numel(intersect(nodes,heavynodes));
    nN=numel(intersect(nodes,normalnodes));
    
    
    if nT==Q %tipEl
        continue
    elseif nH==Q
        heavycount=heavycount+1;
        heavyEl(heavycount)=i;
    elseif nT>0 && nN>0
        neighcount=neighcount+1;
        neighbortipEl(neighcount)=i;
    elseif nT>0 && nH>0 && nN==0
        heavytipcount=heavytipcount+1;
        heavytipEl(heavytipcount)=i;
    else
        normcount=normcount+1;
        normEl(normcount)=i;
    
    end
end


%% Nodvalue
nodvalue=zeros(nnd,1); % 0 = normal
nodvalue(tipnodes)=2;
nodvalue(heavynodes)=1;

if crackedplate==0
    tipEl=[];
    heavyEl=[];
    heavytipEl=[];
    neighbortipEl=[];
    normEl=1:nel;
    tipnodes=[];
    heavynodes=[];
    nodvalue=zeros(nnd,1);
    x1c=0.0001965;
    y1c=0.0001965;
    x2c=0.0001966;
    y2c=0.0001966;
    
end






end






