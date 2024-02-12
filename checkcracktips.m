function checkcracktips(x1c,y1c,x2c,y2c,geom)

s=ismember(geom(:,1),x1c);
t=ismember(geom(:,1),x2c);
u=ismember(geom(:,2),y1c);
v=ismember(geom(:,2),y2c);

h=sum(s+t+u+v);

if h >0
    error('CrackTip Have been located on a Node')
end


end


