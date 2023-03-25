%% Load the cancer data
clear; close all;
load('Pan10.mat');
E = log(1+LUSC_ge);     % BRCA_ge
[num_gene,num_patient] = size(LUSC_ge);


%% TPX2
% find the top 10 related genes to TPX2
seed1 = 'TPX2';
choice = findgene(Gene,seed1);
TPX2_top10 = topgenesmi_fcn(E,Gene,choice);

% iteration loop
TPX2_count = 0; % number of iterations it takes to converge
while 1
    TPX2_count = TPX2_count + 1;

    % Create a fictitious "metagene" whose expression level is the average 
    % of the expression levels of the ten genes that were found
    TPX2_metagene = zeros(1,num_patient);
    for i=1:10
        TPX2_metagene = TPX2_metagene+E(findgene(Gene,TPX2_top10{i}),:);
    end
    TPX2_metagene = TPX2_metagene/10;

    % Find the list of the top ten genes most associated with that metagene 
    F = [E;TPX2_metagene];      % create F by appending the metagene to E
    index_top10 = 2:2+10-1;     % top 10 genes excluding the 1st metagene (ie. 2nd to 11th genes)
    choice = num_gene+1;        % choose the data for metagene (ie. last row of F)
    TPX2_top10_new = topgenesmi_fcn(F,Gene,choice,index_top10);

    % Continue iterating until it converges
    if isequal(TPX2_top10,TPX2_top10_new) == 0
        TPX2_top10 = TPX2_top10_new;
    else
        break;
    end
end

% Compute the mutual information of the final ten genes with the final 
% metagene. The final ten genes are already ranked in order of m.i.
TPX2_top10_new_mi = zeros(1,10);
for i=1:10
    TPX2_top10_new_mi(i) = mi(E(findgene(Gene,TPX2_top10{i}),:),TPX2_metagene);
end

% result: final genes are similar but not exactly the same as before,
% significant discovery that points to the core of expression


%% COL1A2
% find the top 10 related genes to COL1A2
seed2 = 'COL1A2';
choice = findgene(Gene,seed2);
COL_top10 = topgenesmi_fcn(E,Gene,choice);

% iteration loop
COL_count = 0; % number of iterations it takes to converge
while 1
    COL_count = COL_count + 1;

    % Create a fictitious "metagene" whose expression level is the average 
    % of the expression levels of the ten genes that were found
    COL_metagene = zeros(1,num_patient);
    for i=1:10
        COL_metagene = COL_metagene+E(findgene(Gene,COL_top10{i}),:);
    end
    COL_metagene = COL_metagene/10;

    % Find the list of the top ten genes most associated with that metagene 
    F = [E;COL_metagene];       % create F by appending metagene to E
    index_top10 = 2:2+10-1;     % top 10 genes excluding the 1st metagene (ie. 2nd to 11th genes)
    choice = num_gene+1;        % choose the data for metagene (ie. last row of F)
    COL_top10_new = topgenesmi_fcn(F,Gene,choice,index_top10);

    % Continue iterating until it converges
    if isequal(COL_top10,COL_top10_new) == 0
        COL_top10 = COL_top10_new;
    else
        break;
    end
end

% Compute the mutual information of the final ten genes with the final 
% metagene. The final ten genes are already ranked in order of m.i. 
COL_top10_new_mi = zeros(1,10);
for i=1:10
    COL_top10_new_mi(i) = mi(E(findgene(Gene,COL_top10{i}),:),COL_metagene);
end



%% LCP2
% find the top 10 related genes to LCP2
seed3 = 'LCP2';
choice = findgene(Gene,seed3);
LCP2_top10 = topgenesmi_fcn(E,Gene,choice);

% iteration loop
LCP2_count = 0; % number of iterations it takes to converge
while 1
    LCP2_count = LCP2_count + 1;

    % Create a fictitious "metagene" whose expression level is the average 
    % of the expression levels of the ten genes that were found
    LCP2_metagene = zeros(1,num_patient);
    for i=1:10
        LCP2_metagene = LCP2_metagene+E(findgene(Gene,LCP2_top10{i}),:);
    end
    LCP2_metagene = LCP2_metagene/10;

    % Find the list of the top ten genes most associated with that metagene 
    F = [E;LCP2_metagene];      % create F by appending the metagene info to E
    index_top10 = 2:2+10-1;     % top 10 genes excluding the 1st metagene (ie. 2nd to 11th genes)
    choice = num_gene+1;        % choose the data for metagene (ie. last row of F)
    LCP2_top10_new = topgenesmi_fcn(F,Gene,choice,index_top10);

    % Continue iterating until it converges
    if isequal(LCP2_top10,LCP2_top10_new) == 0
        LCP2_top10 = LCP2_top10_new;
    else
        break;
    end
end

% Compute the mutual information of the final ten genes with the final 
% metagene. The final ten genes are already ranked in order of m.i.
LCP2_top10_new_mi = zeros(1,10);
for i=1:10
    LCP2_top10_new_mi(i) = mi(E(findgene(Gene,LCP2_top10{i}),:),LCP2_metagene);
end


%% Scatter plot of 3 genes for the max-strength signature
% LCP2 has the highest strength = 0.8902 (lowest/tenth of m.i.)
LCP2_g1 = LCP2_top10_new_sort(1);LCP2_g1_loc = find(strcmp(Gene,LCP2_g1) == 1);
LCP2_g2 = LCP2_top10_new_sort(2);LCP2_g2_loc = find(strcmp(Gene,LCP2_g2) == 1);
LCP2_g3 = LCP2_top10_new_sort(3);LCP2_g3_loc = find(strcmp(Gene,LCP2_g3) == 1);
 
sc3(LCP2_g1_loc, LCP2_g2_loc, LCP2_g3_loc, E);
xlabel(LCP2_g1);ylabel(LCP2_g2);


%% Scatter plot for other signatures
% COL1A2
COL_g1 = COL_top10_new_sort(1);COL_g1_loc = find(strcmp(Gene,COL_g1) == 1);
COL_g2 = COL_top10_new_sort(2);COL_g2_loc = find(strcmp(Gene,COL_g2) == 1);
COL_g3 = COL_top10_new_sort(3);COL_g3_loc = find(strcmp(Gene,COL_g3) == 1);
sc3(COL_g1_loc, COL_g2_loc, COL_g3_loc, E);
xlabel(COL_g1);ylabel(COL_g2);

% TPX2
TPX2_g1 = TPX2_top10_new_sort(1);TPX2_g1_loc = find(strcmp(Gene,TPX2_g1) == 1);
TPX2_g2 = TPX2_top10_new_sort(2);TPX2_g2_loc = find(strcmp(Gene,TPX2_g2) == 1);
TPX2_g3 = TPX2_top10_new_sort(3);TPX2_g3_loc = find(strcmp(Gene,TPX2_g3) == 1);
sc3(TPX2_g1_loc, TPX2_g2_loc, TPX2_g3_loc, E);
xlabel(TPX2_g1);ylabel(TPX2_g2);


%% Survival analysis of breast cancer data
load('TCGA_PanCan12.mat');
E = log(1+BRCA_ge); 
Time = BRCA_time; 
Status = BRCA_status;
[num_gene,num_patient] = size(BRCA_ge);

% Compute concordance index (CI) for all genes
CI = zeros(1,num_gene);
for i=1:num_gene
    CI(i) = concordanceindex(E(i,:),Time,Status);
end

% Rank all genes in terms of their concordance index (CI)
[CI_sort,I] = sort(CI);
gene_sort = Gene(I);
top5_protective_genes = gene_sort(1:5);


%% KM curve for LOC100128977
most_protective_loc = findgene(Gene,'LOC100128977');
xx = find(E(most_protective_loc,:) > median(E(most_protective_loc,:)));
yy = find(E(most_protective_loc,:) <= median(E(most_protective_loc,:)));

figure;
kaplan_meier(Time(xx),Status(xx),'b');  % high LOC100128977 expression
hold on;
kaplan_meier(Time(yy),Status(yy),'r');  % low LOC100128977 expression
title('Kaplan-Meier Plot for LOC100128977');
xlabel('Time (days)');
ylabel('Percent survival of patients (%)');
legend('subpopulation with high LOC100128977 expression','subpopulation with low LOC100128977 expression');


%% Heat map
E = log(1 + LUSC_ge); 

% Initiate signature genes' location
sig1 = zeros(1,10);
sig2 = sig1;
sig3 = sig1;

% Location of the genes of each signature
for i = 1:10
    sig1(i) = findgene(Gene,TPX2_top10_new_sort{i});
    sig2(i) = findgene(Gene,COL_top10_new_sort{i});
    sig3(i) = findgene(Gene,LCP2_top10_new_sort{i});
end

% function "mean" outputs the mean values of each column 
[~,I1] = sort(mean(E(sig1,:)),'descend');   
[~,I2] = sort(mean(E(sig2,:)),'descend');   
[~,I3] = sort(mean(E(sig3,:)),'descend');   

selected_genes=[sig1 sig2 sig3];
selected_samples=[I1(1:10) I2(1:10) I3(1:10)];
G = E(selected_genes,selected_samples);
clustergram(G,'Standardize','row','RowLabels',Gene(selected_genes),'ColumnLabels',selected_samples);


