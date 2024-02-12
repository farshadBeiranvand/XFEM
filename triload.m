function [nodal_lodes]=triload(Q,p0,pL,len,Nelement)
%(p0= p at x,y=0) (pL= p at x,y=L,w)
b=len/Nelement;

if Q==4
error('the distributed load (load factors) didnt defind for Q=4')
end

mp=(pL-p0)/len;  %f= mp*x+p0


elements_load=zeros(Nelement,3);

for i=1:Nelement
	pbeg=mp*((i-1)*b)+p0;  %pbegining
	pend=mp*(i*b)+p0;      %pend

		if mp>0

			if pbeg+pend>0  %compressive force
			p_uniformload=[1/6*pbeg 2/3*pbeg 1/6*pbeg];
			ptri=pend-pbeg;
			p_triload=[0 1/3*ptri 1/6*ptri];
			elements_load(i,:)=p_uniformload+p_triload;
			elseif pbeg+pend<0 %traction
			p_uniformload=[1/6*pend 2/3*pend 1/6*pend];
			ptri=pbeg-pend;
			p_triload=[1/6*ptri 1/3*ptri 0];
			elements_load(i,:)=p_uniformload+p_triload;
			end


		elseif mp<0

			if pbeg+pend>0 
			p_uniformload=[1/6*pend 2/3*pend 1/6*pend];
			ptri=pbeg-pend;
			p_triload=[1/6*ptri 1/3*ptri 0];
			elements_load(i,:)=p_uniformload+p_triload;
			elseif pbeg+pend<0 
			p_uniformload=[1/6*pbeg 2/3*pbeg 1/6*pbeg];
			ptri=pend-pbeg;
			p_triload=[0 1/3*ptri 1/6*ptri];
			elements_load(i,:)=p_uniformload+p_triload;
			end

		end
end %for


Eload=reshape(elements_load.',1,[]);
Nnodes=2*Nelement+1;
nodal_lodes=zeros(Nnodes,1);

nodal_lodes(1)=Eload(1);
nodal_lodes(end)=Eload(end);

q=2;
for k=2:Nnodes-1

	if rem (k,2)==0
	nodal_lodes(k)=Eload(q);
	q=q+1;
	elseif rem (k,2)==1
	nodal_lodes(k)=Eload(q)+Eload(q+1);
	q=q+2;
	end
end




end
	







