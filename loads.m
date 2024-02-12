
function [force]=loads(Q,gdofIN,geom,w,L,NYE,NXE,Nxxlo,Nxylo,Nyylo,Nxxnl,Nxynl,Nyynl)
%force vector
%compressive force
geom=vpa(geom);
force=zeros(gdofIN,1);
left=find(geom(:,1)==0);
right=find(geom(:,1)==2*w);
down=find(geom(:,2)==0);
up=find(geom(:,2)==2*L);

px=Nxxlo*2*L/NYE;
py=Nyylo*2*w/NXE;

u_left=2*left-1;
u_right=2*right-1;
v_down=2*down;
v_up=2*up;

if Q==8
force(u_left(1:2:end))=px/3;
force(u_left(2:2:end))=px*2/3;
force(u_left(1))=px/6;
force(u_left(end))=px/6;

force(u_right(1:2:end))=-px/3;
force(u_right(2:2:end))=-px*2/3;
force(u_right(1))=-px/6;
force(u_right(end))=-px/6;

force(v_down(1:2:end))=py/3;
force(v_down(2:2:end))=py*2/3;
force(v_down(1))=py/6;
force(v_down(end))=py/6;

force(v_up(1:2:end))=-py/3;
force(v_up(2:2:end))=-py*2/3;
force(v_up(1))=-py/6;
force(v_up(end))=-py/6;



elseif Q==4
force(u_left)=px;
force(u_left(1))=px/2;
force(u_left(end))=px/2;

force(u_right)=-px;
force(u_right(1))=-px/2;
force(u_right(end))=-px/2;

force(v_down)=py;
force(v_down(1))=py/2;
force(v_down(end))=py/2;

force(v_up)=-py;
force(v_up(1))=-py/2;
force(v_up(end))=-py/2;

end
end



