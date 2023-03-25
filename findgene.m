function index=findgene(genes,s)
for i=1:length(genes)
    if isequal(genes{i},s)
        index=i;
    end
end
end