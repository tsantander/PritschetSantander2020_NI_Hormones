function [correlation] = bgt_correlationMatrix(filteredSeries, sigTest, fig)
% Compute functional association via Pearson product-moment correlations.
%
% FORMAT [correlation] = bgt_correlationMatrix(filteredSeries, sigTest, fig)
%
% REQUIRED INPUT:
%   filteredSeries
%       m x n matrix of temporally-filtered data, where m = number of
%       timepoints and n = number of regions of interest (ROIs) in the 
%       network.
%
%   sigTest
%       Indicates a method for calculating p-values (i.e. the probability
%       of obtaining a given correlation value, assuming the null of no
%       association between brain regions). Broadly, these fall into
%       parametric or surrogate (nonparametric) approaches. 'Parametric'
%       estimation is fastest. Nonparametric estimation (i.e. permutation 
%       testing or surrogate amplitude-adjusted Fourier transforms) allows
%       us to obtain empirical null distributions, which I generally
%       prefer, but it comes at computational cost. Specify (in single quot
%       quotes):
%           'parametric'        -   Standard estimation of p-values for
%                                   correlation, based on a transformation
%                                   to t-statistics.
%           'NPT'               -   Uses random permutations of the
%                                   timeseries to build empirical null
%                                   distributions. This is slightly less
%                                   robust than AAFT/IAAFT.
%           'AAFT'              -   Uses amplitude-adjusted Fourier
%                                   transforms to generate surrogate data
%                                   for each pairwise comparison. In brief,
%                                   AAFT produces phase-shuffled null
%                                   samples while preserving the original
%                                   amplitude distribution of the data.
%                                   Note, however, that these
%                                   transformations may alter the assumed
%                                   linear structure of the signal..
%           'IAAFT'             -   Uses iterative AAFT to generate
%                                   surrogate data. This will take longer
%                                   than standard AAFT but provides a
%                                   'better' approximation of the original
%                                   signal's autocorrelation function in
%                                   addition to its amplitude distribution.
%
%   fig
%       Indicates whether or not to display the correlation heatmap. Enter
%       1 for YES or 0 for NO.
%
% OUTPUT:
%   correlation
%       Structure array with the following fields:
%           .uncorrected
%               .associationMatrix     -    Symmetric n x n matrix whose
%                                           elements (i,j) indicate the
%                                           correlation between the ith
%                                           and jth ROIs in your network.
%               .pValues               -    Symmetric n x n matrix whose
%                                           elements (i,j) indicate the
%                                           probability of obtaining the
%                                           correlation value contained
%                                           in associationMatrix(i,j) by
%                                           chance (if the null of no
%                                           association is true).
%__________________________________________________________________________
%
% This function will compute the extent of co-activity between pairs of
% brain regions using Pearson correlation.
%__________________________________________________________________________
%
% BRAIN GRAPHS: A toolbox for graph theoretic analyses of fMRI data, v1.03
% Author:
%   Tyler Santander (t.santander@psych.ucsb.edu)
%   Institute for Collaborative Biotechnologies
%   Department of Psychological & Brain Sciences
%   University of California, Santa Barbara
%   December 2018
%__________________________________________________________________________
          
% Compute Pearson correlation between all nodal timeseries.
%--------------------------------------------------------------------------

    nROI        = size(filteredSeries,2);
    nEdge       = (nROI^2 - nROI)/2;
    correlation = [];
        
    startCor = tic;
            
    disp(['|| Estimating correlation between ' num2str(nEdge) ' network edges...']);
    
    switch sigTest
        
        case 'parametric'
    
            [associationMatrix, pValues] = corrcoef(filteredSeries);
            
            % Clear diagonals.
            
                associationMatrix(1:nROI+1:end) = 0;
                pValues(1:nROI+1:end)           = 0;
            
        case {'NPT', 'AAFT', 'IAAFT'}
            
            associationMatrix = zeros(nROI);
            pValues           = zeros(nROI);
            
            for iROI = 1:nROI
                
                for jROI = 1:nROI
                    
                    if (iROI <= jROI)   % Skip diagonal and lower triangle.
                        
                        continue
                        
                    else                % Otherwise, estimate correlation.
                    
                        [rTrue, ~] = corrcoef(filteredSeries(:,iROI), filteredSeries(:,jROI));
                    
                        associationMatrix(iROI,jROI) = rTrue(1,2);
                        
                        % Generate 1000 shuffled samples of a
                        % timeseries (either via permutation or (I)AAFT.
                                
                        nIter   = 1000;
                        randCor = zeros(nIter,1);
                        
                        switch sigTest
                            
                            case 'NPT'
                                
                                temp  = filteredSeries(:,jROI);
                                yPerm = zeros(tsLength, nIter);
                                
                                for iPerm = 1:nIter
                                    yPerm(:,iPerm) = temp(randperm(tsLength));
                                end
                                
                                clear temp
                                
                                % Loop over permutations.
                                
                                for iPerm = 1:nIter
                                    
                                    [rPerm, ~]     = corrcoef(filteredSeries(:,iROI), yPerm(:,iPerm));
                                    randCor(iPerm) = rPerm(1,2);
                                    
                                end
                                
                            case 'AAFT'
                                
                                xPerm = AAFT(filteredSeries(:,iROI), nIter);
                                yPerm = AAFT(filteredSeries(:,jROI), nIter);
                                
                                % Loop over surrogates.
                                
                                for iPerm = 1:nIter
                                    
                                    [rPerm, ~]     = corrcoef(xPerm(:,iPerm), yPerm(:,iPerm));
                                    randCor(iPerm) = rPerm(1,2);
                                
                                end
                                
                                clear xPerm yPerm
                            
                            case 'IAAFT'
                                
                                xPerm = IAAFT(filteredSeries(1:tsLength,iROI), nIter);
                                yPerm = IAAFT(filteredSeries(1:tsLength,jROI), nIter);
                                
                                % Loop over surrogates.
                                
                                for iPerm = 1:nIter
                                    
                                    [rPerm, ~]     = corrcoef(xPerm(:,iPerm), yPerm(:,iPerm));
                                    randCor(iPerm) = rPerm(1,2);
                                
                                end
                                
                                clear xPerm yPerm
                                
                        end
                        
                        % Get an empirical p-value based on the proportion
                        % of permuted correlations that were greater than 
                        % or equal to the 'true' estimate for this pair of
                        % timeseries (two-tailed).
                        
                            pValues(iROI,jROI) = sum(double(abs(randCor) >= abs(rTrue)))/nIter;
                    
                    end
                    
                end
                
                associationMatrix = associationMatrix + associationMatrix.';
                pValues           = pValues + pValues.';
                
            end
            
    end
            
    endCor = toc(startCor);
    disp(['|| Correlation between ' num2str(nEdge) ' network edges computed in ' num2str(endCor/60) ' minutes']);
            
    correlation.uncorrected.associationMatrix = associationMatrix;
    correlation.uncorrected.pValues           = pValues;
            

    % Display correlation heatmap (optional, set 'fig' argument to 1 if
    % YES).
    
        if fig
            
            figure; image(correlation.uncorrected.associationMatrix .* 64);
            
        end

        
% Save correlation structure.
%--------------------------------------------------------------------------

    save('correlation', 'correlation');

end