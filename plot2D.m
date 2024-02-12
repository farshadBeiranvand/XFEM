function plot2D(NXE,NYE,heavyEl,neighbortipEl,heavytipEl,tipEl,x1c,x2c,y1c,y2c,connec,geom)
figure('Color','white')

enrichvalue=ones(NXE*NYE,1); %value of enrichment for each element   1:normal 2:heavy 3:neighbortip 4:heavytip 5:tip



% #plot plate with crack and enrichment elements

facevertexdata=zeros(NXE*NYE,3);
for ii=1:NXE*NYE
switch enrichvalue(ii,1)
case 1
facevertexdata(ii,:)=[1 1 1]; %white for normal
case 2
facevertexdata(ii,:)=[1 1 0]; %yellow for heavy 
case 3
facevertexdata(ii,:)=[.8 1 1]; %Turquoise for neighbortip
case 4
facevertexdata(ii,:)=[0 1 .2]; %green for heavytip
case 5
facevertexdata(ii,:)=[0 0 1]; %blue for tip
end
end
patch('Faces', connec, 'Vertices', geom, 'FaceVertexCData',facevertexdata,'Facecolor','flat','Marker','o','MarkerSize',3);



end
