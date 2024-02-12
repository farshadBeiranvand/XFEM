function newgeom=newgeom(displacements,nodvalue,heavynodes,tipnodes,geom,nnd);
M=5e5; % magnify crack displacements
disp=displacements;
[a,b]=size(geom);
newgeom=geom(a,b);
nHG=numel(heavynodes);
for i=1:nnd
    m=nodvalue(i);
    if m==0
        newgeom(i,1)=geom(i,1)+disp(2*i-1);
        newgeom(i,2)=geom(i,2)+disp(2*i);
    elseif m==1
        F=find(heavynodes==i);
        newgeom(i,1)=geom(i,1)+disp(2*i-1)+disp(2*nnd+2*F-1)*M;
        newgeom(i,2)=geom(i,2)+disp(2*i)+disp(2*nnd+2*F)*M;
    elseif m==2
        T=find(tipnodes==i);
        ddf=2*nnd+2*nHG+8*T;
        newgeom(i,1)=geom(i,1)+disp(2*i-1)+M/10*(disp(ddf-7)+disp(ddf-5)+disp(ddf-3)+disp(ddf-1));
        newgeom(i,2)=geom(i,2)+disp(2*i)+M/10*(disp(ddf-6)+disp(ddf-4)+disp(ddf-2)+disp(ddf));
    end
end


end