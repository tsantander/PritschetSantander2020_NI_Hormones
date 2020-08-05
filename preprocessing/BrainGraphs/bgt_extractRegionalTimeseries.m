function [timeSeries] = bgt_extractRegionalTimeseries(data4D, ROIs, summaryMeasure)
% Regional timeseries extraction.
%
% FORMAT [timeSeries] = bgt_extractRegionalTimeseries(data4D, ROIs, summaryMeasure)
%
% REQUIRED INPUT:
%   data4D
%       String specifying file/path information for 4D (x,y,z + time)
%       image volume.
%
%   ROIs
%       String specifying file/path information for 3D image volume (e.g.
%       a cortical parcellation), where each voxel contains a scalar label 
%       indicating ROI membership. NOTE: Assumes this image has been
%       coregistered to the volumes specified in data4D (so they are in the
%       same space / have the same dimensions).
%
%   summaryMeasure
%       Summary statistic of choice. Specify (in single quotes) either:
%           'mean'      - Extract mean signal within each ROI.
%           'eigen1'    - Extract first eigenvariate within each ROI. 
%                         Especially appropriate if dealing with large
%                         anatomical ROIs, as the voxelwise response may
%                         be highly heterogeneous. The first eigenvariate
%                         will capture the dimension along which maximal
%                         variance is accounted for (so it's effectively
%                         a 'weighted' mean).
%
% OUTPUT:
%   timeSeries
%       m x n matrix where m = number of timepoints and n = number of nodes
%       (i.e. ROIS) in the network.
%__________________________________________________________________________
%
% This function will extract regional timeseries from a 4D image volume.
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

% Load image volumes into memory.
%--------------------------------------------------------------------------

    % 4D functional images.
    
        V         = load_nii(data4D);
        data4D    = V.img;
        [x,y,z,t] = size(data4D);
        data2D    = double(reshape(data4D, x*y*z, t)');
    
        clear V data4D
    
    % 3D ROI image.
    
        M         = load_nii(ROIs);
        mask3D    = M.img;
        [x,y,z]   = size(mask3D);
        mask2D    = reshape(mask3D, 1, x*y*z);
    
        clear M mask3D
    
    % Quick check to make sure images are same size.
    
        if (size(data2D,2) ~= length(mask2D))
            error('- ERROR: Image volumes have unequal spatial dimensions! -');
        end
    
        
% Initialize timeSeries matrix.   
%--------------------------------------------------------------------------

    timeSeries = zeros(t, max(mask2D));
    
    
% Begin looping over ROIs, extract regional signal as specified in
% summaryMeasure.
%--------------------------------------------------------------------------

    disp('|| Extracting regional timeseries...');
    
    startExtract = tic;

    for iROI = 1:max(mask2D)
        
        dataROI = data2D(:, mask2D == iROI);
        
        switch summaryMeasure
            
            case 'mean'
                
                % Summarize regional response in terms of mean.
                
                timeSeries(:,iROI) = mean(dataROI,2);
                
            case 'eigen1'   
                
                % Compute regional response in terms of first eigenvariate.
                % Code taken from spm_regions, with thanks to Karl Friston.
                
                [m,n] = size(dataROI);
                
                if m > n
                    
                    [v, s, v] = svd(dataROI'*dataROI);
                    s         = diag(s);
                    v         = v(:,1);
                    u         = dataROI*v/sqrt(s(1));
                    
                else
                    
                    [u, s, u] = svd(dataROI*dataROI');
                    s         = diag(s);
                    u         = u(:,1);
                    v         = dataROI'*u/sqrt(s(1));
                    
                end
                
                d                  = sign(sum(v));
                u                  = u*d;
                v                  = v*d;
                timeSeries(:,iROI) = u*sqrt(s(1)/n);
                
            otherwise
                
                error('- ERROR: Invalid summary function specified! -');
        end
        
    end
    
    endExtract = toc(startExtract);
    
    disp(['|| Timeseries extracted from ' num2str(max(mask2D)) ' ROIs in ' num2str(endExtract) ' seconds']);
    
    
% Save timeseries.
%--------------------------------------------------------------------------

    save('timeSeries','timeSeries');
    
end
        
        
        