function padded = padmap(map, pad)
padded=zeros(size(map,1)+(size(map,1)-1)*pad , size(map,2));
for i=1:size(map,2)-1
    st = map(i+1,:);
    go = map(i,:);
    phi = linspace(0,1,pad+2);
    ilo = (i-1)*(pad+2)+1;
    ihi = i*(pad+2);
    padded(ilo:ihi,:)=(st-go).*phi' + go;
end

end