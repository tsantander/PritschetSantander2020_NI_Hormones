function bgt_globalNorm(timeSeries, tendencyMeasure)
% Global intensity normalization.
%
% FORMAT bgt_globalNorm(timeSeries, tendencyMeasure)
%
% REQUIRED INPUT:
%   timeSeries
%       Data can be specified here in one of two forms:
%           1) String specifying file/path information for 4D
%              (x,y,z + time) image volume.
%              NOTE: must be a character array, not a cell array.
%           2) m x n matrix where m = number of TRs and n = number of
%              voxels or regions of interest (ROIs) in your data.
%
%   tendencyMeasure
%       Central tendency measure of choice. Specify (in single quotes)
%       either:
%           'mean'      - Scale by grand mean.
%           'median'    - Scale by grand median.
%
% OUTPUT:
%   No explicit output. If timeSeries was entered as an image volume, new
%   images will be written to the working directory; if timeSeries was
%   entered as a matrix, a new .mat file will be written to the working
%   directory.
%__________________________________________________________________________
%
% This function will scale the voxelwise/ROIwise timeseries according to
% either the grand mean or grand median. Global intensity normalization
% may reduce the risk of mischaracterizing functional associations
% between brain regions as a result of transient fluctuations in signal
% intensity across space/time. Also helps for subsequent cross-subject
% comparison by ensuring all timeseries have same grand mean/median.
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

% Check timeseries input argument. If matrix, rename. Otherwise extract
% voxelwise timeseries from preprocessed EPIs.
%--------------------------------------------------------------------------

    tsType = ischar(timeSeries);

    switch tsType
        
        case 0
            
            if (iscell(timeSeries))
                
                error('- ERROR: timeSeries argument must be a character array or matrix! -');
                
            else
            
                tsMatrix = timeSeries;
                
            end
    
        case 1
        
            % Load image volumes into memory.
            
                V               = load_nii(timeSeries);
                data4D          = V.img;
                dataHeader      = V.hdr;
                [path,name,ext] = fileparts(timeSeries);  
    
            % Compress out non-brain voxels (implicit mask) for 
            % computational efficiency later.
    
                [x,y,z,t]      = size(data4D);
                data2D         = reshape(data4D, x*y*z, t);
                mInd           = all(data2D == 0, 2);
                data2D(mInd,:) = [];
                tsMatrix       = double(data2D)';
                
                clear V data2D data4D
                
        otherwise
            
            error('- ERROR: timeSeries argument must be a character array or matrix! -');
    
    end
    
    
% Scale the data across space and time.
%--------------------------------------------------------------------------

    disp(['|| Scaling data according to the grand ' tendencyMeasure '...']);

    switch tendencyMeasure
        
        case 'mean'
        
            globalValue = mean(tsMatrix(:));
    
        case 'median'
        
            globalValue = median(tsMatrix(:));
            
        otherwise
            
            error('- ERROR: Invalid tendency measure specified! -');

    end
    
    scaleFactor    = 1000/globalValue;
    tsMatrixScaled = tsMatrix*scaleFactor;
    
    disp('|| Finished scaling');
    
    
% Write output.
%--------------------------------------------------------------------------

    switch tsType
        
        case 0
            
            save('globalNormSeries','tsMatrixScaled');
            
        case 1
                        
            data2D              = zeros(t, x*y*z);
            data2D(:,mInd == 0) = tsMatrixScaled;
            data4D              = reshape(data2D', [x y z t]);
            
            nii.img             = data4D;
            nii.hdr             = dataHeader;
            
            save_nii(nii, [path '/g' name ext]);
            
    end
    
end