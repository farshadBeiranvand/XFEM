%remove big displacement
function Wnodes=removeD(Wnodes)
b=sort(Wnodes);
c=unique(b);
m1=c(end);
m2=c(end-1);
m3=c(end-2);
m4=c(end-3);
for i=max(size(Wnodes))
    k=Wnodes(i);
    if (k==m1|| k==m2|| k==m3)
     Wnodes(i)=m4;
    end
end