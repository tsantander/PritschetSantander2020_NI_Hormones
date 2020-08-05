function [P,thresh_val] = mod_allegiance(C,toThresh)
% given a set of partitions, this function returns module allegiance matrix
% of probabilities that any two nodes will be associated together in the
% same module (community) over the given partitions
% INPUTS: C, a pxn matrix containing p partitions of an n-node network
%         toThresh, a boolean value indicating whether to threshold output
% OUTPUTS: P, an nxn matrix where entry Pij contains the fraction of input
%           partitions in which nodes i & j were grouped in the same module
%          thresh_val, the scalar value below which P was thresholded
%                (if toThresh==0, thresh_val = min value of P)

if nargin<2
    toThresh = 0;
end

[p,n] = size(C);

P = zeros(n);
if toThresh
    Prand = zeros(n);
   
    % create random partitions for significance threshold
    Crand = zeros(size(C));
    for i=1:p
        pr = randperm(n);      % create random partitions
        Crand(i,:) = C(i,pr);  % Crand, p random assignments of n nodes
    end

end

% create nxn count of partitions in which nodepairs are in the same module
    for ii = 1:n          
        for jj = ii:n;
            
            % actual version
            P(ii,jj) = sum(C(:,ii)==C(:,jj));
            
            % random versions
            if toThresh
                Prand(ii,jj) = sum(Crand(:,ii)==Crand(:,jj));
            end
            
        end % end loop over nodepairs
    end % end loop over second of nodepairs    

% symmetrise P (no need to symmetrise Prand; we only need its max value)
P = P + triu(P,1)';

% threshold, if desired
if toThresh
    thresh_val = max(max(triu(Prand,1)));
    P(P<=thresh_val) = 0;
else
    thresh_val = min(min(P));
end

% normalize P by number of partitions
P = P./p;

end