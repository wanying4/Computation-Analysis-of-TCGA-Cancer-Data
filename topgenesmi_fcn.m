function [ top10,top10_mi ] = topgenesmi_fcn(E,Gene,choice,index_top10)
%TOPGENESMI_FCN finds the top ten genes in terms of their mutual information 
%with the gene defined by Gene(choice).
%   Inputs:
%       E is the expression matrix, rows are genes, cols are patient 
%       biopsy's result.
%       Gene is an array that contains a list of gene names.
%       choice is the index of the gene in interest in the array Gene.
%       index_top10 by default should be 1:10.

[N,~] = size(E);  % [N,M]=size(E);
vect = -ones(N,1);
y = E(choice,:)';

for i=1:N
    x = E(i,:);
    if sum(x)>0 %exclude genes with zeros throughout
        vect(i)= mi(x,y);
    end
    if rem(i,1000)==0
        i;
    end
end

% sort vect
[Y,I]=sort(vect,'descend');
top10_mi = Y(2:11);

% If the index_top10 variable exist (ie. defined) then use that index for
% top 10 genes. Else the default will be the 1st to 10th gene.
if exist('index_top10','var') == 1
    top10 = Gene(I(index_top10));
else
    top10 = Gene(I(1:10));
end

end