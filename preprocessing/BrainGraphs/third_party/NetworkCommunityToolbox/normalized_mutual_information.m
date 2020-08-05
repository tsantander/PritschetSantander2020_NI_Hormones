function [nmi]=normalized_mutual_information(partition1,partition2,varargin)

% NORMALIZED_MUTUAL_INFORMATION - This function computes the normalized mutual
%         information (normalization taken from from Danon 2005, note other normalizations
%         exist) which is a measure of similarity between two different partitions of 
%         community structure.  The measure is bounded in [0 1] where 1 implies perfect
%         agreement between partitions.  In cases where the number of nodes
%         divided by the number of clusters <~ 100, the adjusted value should be
%         used to correct for chance (see Vinh et al 2010).
% 
% Inputs: 
%         partition1 - vector containing community assignment for each node in first partition
% 
%         partition2 - vector containing community assignment for each node in second partition
% 
%         varargin - set equal to 'adjusted' if want to correct for chance and compute the AMI 
%               (use this option if the number of nodes diveded by the number of communities <~ 100.  
%               Default is the unadjusted NMI (set varaible equal to 'unadjusted' or leave empty).
% 
% Outputs: 
%         nmi - normalized mutual information.  Gives measure of
%               similarity between partitons.  can be adjusted for chance 
%               (AMI- adjusted normalized mutual information - see Vihn 2010)
% 
% Written by Sarah Feldt Muldoon

% References: Danon L, Diaz-Guilera A, Duch J, Arenas A (2005) Comparing
%         community structure identification. J Stat Mech:P09008.; Vinh NX, Epps J, 
%         Bailey J (2010) Information Theoretic Measures for Clusterings Comparison: 
%         Variants, Properties, Normalization and Correction for Chance. The Journal
%         of Machine Learning Research.


%if no keyword, calculate the unadjusted NMI
if isempty(varargin)
    varargin{1} = 'unadjusted';
end

num_nodes = length(partition1);
comm1_array = unique(partition1);
comm2_array = unique(partition2);
num_comm1 = length(comm1_array);
num_comm2 = length(comm2_array);
contingency_matrix = zeros(num_comm1,num_comm2);

%make contingency table (also called confusion_matrix)
for i = 1:num_comm1
    comm1 = comm1_array(i);
    for j = 1:num_comm2
        comm2 = comm2_array(j);
        contingency_matrix(i, j) = sum(partition1 == comm1 & partition2 == comm2);
    end
end

sum_over_i = sum(contingency_matrix,2);
sum_over_j = sum(contingency_matrix,1);
total_sum = sum(sum(contingency_matrix));

%compute mutual information
mi_matrix = contingency_matrix .* log((contingency_matrix .* total_sum) ./ (sum_over_i * sum_over_j));
mi = sum(mi_matrix(isfinite(mi_matrix))) / total_sum;

%compute normalization terms (we normalize by (h1+h2)/2 as in Danon 2005)
h1 = sum_over_i .* log(sum_over_i ./ total_sum);
h1 = -sum(h1(isfinite(h1))) / total_sum;
h2 = sum_over_j .* log(sum_over_j ./ total_sum);
h2 = -sum(h2(isfinite(h2))) / total_sum;

if strcmp(varargin{1}, 'adjusted') == 1
    %compute expected mutual information (Vihn 2010)
    expected_mi = 0;
    for i = 1:num_comm1
        for j=1:num_comm2
            k_min = max([sum_over_i(i) + sum_over_j(j) - total_sum, 1]);
            k_max = min([sum_over_i(i), sum_over_j(j)]);
            for k = k_min:k_max
                term1 = k / total_sum;
                term2 = log(total_sum * k / (sum_over_i(i) * sum_over_j(j)));
                term3_numerator = factorial(sum_over_i(i)) * factorial(sum_over_j(j)) * factorial(total_sum - sum_over_i(i)) * factorial(total_sum - sum_over_j(j));
                term3_denominator = factorial(total_sum) * factorial(k) * factorial(sum_over_i(i) - k) * factorial(sum_over_j(j) - k) * factorial(total_sum - sum_over_i(i) - sum_over_j(j) + k);
                emi_term = term1 * term2 * term3_numerator / term3_denominator;
                if isfinite(emi_term)
                    expected_mi = expected_mi + emi_term;
                end
            end
        end
    end
    %compute AMI_sum as in Vihn 2010
    nmi = (mi - expected_mi) / (.5 * (h1 + h2) - expected_mi);
elseif strcmp(varargin{1}, 'unadjusted') == 1
    %compute NMI as in Danon 2005 (NMI_sum in Vihn 2010)
    nmi = 2 * mi / (h1 + h2);
end



