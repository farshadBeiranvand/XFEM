
function [force]=loads2(Q,gdofIN,geom,w,L,NYE,NXE,Nxxlo,Nxylo,Nyylo,Lfx,Lfy)



geom=vpa(geom);
force=zeros(gdofIN,1);
left=find(geom(:,1)==0);
right=find(geom(:,1)==2*w);
down=find(geom(:,2)==0);
up=find(geom(:,2)==2*L);

u_left=2*left-1;
u_right=2*right-1;
v_down=2*down;
v_up=2*up;

px=Nxxlo*2*L/NYE;
switch Lfx
    case 0
        if Q==8
            force(u_left(1:2:end))=px/3;
            force(u_left(2:2:end))=px*2/3;
            force(u_left(1))=px/6;
            force(u_left(end))=px/6;
            
            force(u_right(1:2:end))=-px/3;
            force(u_right(2:2:end))=-px*2/3;
            force(u_right(1))=-px/6;
            force(u_right(end))=-px/6;
            
        elseif Q==4
            force(u_left)=px;
            force(u_left(1))=px/2;
            force(u_left(end))=px/2;
            
            force(u_right)=-px;
            force(u_right(1))=-px/2;
            force(u_right(end))=-px/2;
        end
        
    case 0.5
        [nodal_lodes]=triload(Q,px,0.5*px,L,NYE);
        force(u_left)=nodal_lodes;
        force(u_right)=-nodal_lodes;
        
    case 1
        [nodal_lodes]=triload(Q,px,0,L,NYE);
        force(u_left)=nodal_lodes;
        force(u_right)=-nodal_lodes;
        
    case 1.5
        n1=fix(.75*NYE);
        n2=NYE-n1;
        [L1]=triload(Q,px,0,L*n1/NYE,n1);
        [L2]=triload(Q,0,-0.5*px,L*n2/NYE,n2);
        L1(end)=[];
        nodal_lodes=vertcat(L1,L2);
        force(u_left)=nodal_lodes;
        force(u_right)=-nodal_lodes;
        
    case 2
        n1=fix(.5*NYE);
        n2=NYE-n1-1;
        [L1]=triload(Q,px,0,L*n1/NYE,n1);
        [L2]=triload(Q,0,-px,L*n2/NYE,n2);
        L1=[L1;0];
        nodal_lodes=vertcat(L1,L2);
        force(u_left)=nodal_lodes;
        force(u_right)=-nodal_lodes;
end

py=Nyylo*2*w/NXE;
switch Lfy
    case 0
        
        
        if Q==8
            force(v_down(1:2:end))=py/3;
            force(v_down(2:2:end))=py*2/3;
            force(v_down(1))=py/6;
            force(v_down(end))=py/6;
            
            force(v_up(1:2:end))=-py/3;
            force(v_up(2:2:end))=-py*2/3;
            force(v_up(1))=-py/6;
            force(v_up(end))=-py/6;
            
        elseif Q==4
            force(v_down)=py;
            force(v_down(1))=py/2;
            force(v_down(end))=py/2;
            
            force(v_up)=-py;
            force(v_up(1))=-py/2;
            force(v_up(end))=-py/2;
        end
        
    case 0.5
        [nodal_lodes]=triload(Q,py,0.5*py,w,NXE);
        force(v_down)=nodal_lodes;
        force(v_up)=-nodal_lodes;
        
    case 1
        [nodal_lodes]=triload(Q,py,0,w,NXE);
        force(v_down)=nodal_lodes;
        force(v_up)=-nodal_lodes;
        
    case 1.5
        n1=fix(.75*NXE);
        n2=NXE-n1;
        [L1]=triload(Q,py,0,w*n1/NXE,n1);
        [L2]=triload(Q,0,-0.5*py,w*n2/NXE,n2);
        L1(end)=[];
        nodal_lodes=vertcat(L1,L2);
        force(v_down)=nodal_lodes;
        force(v_up)=-nodal_lodes;
        
    case 2
        n1=fix(.5*NXE);
        n2=NXE-n1-1;
        [L1]=triload(Q,py,0,w*n1/NXE,n1);
        [L2]=triload(Q,0,-py,w*n2/NXE,n2);
        L1=[L1;0];
        nodal_lodes=vertcat(L1,L2);
        force(v_down)=nodal_lodes;
        force(v_up)=-nodal_lodes;
        
end



end