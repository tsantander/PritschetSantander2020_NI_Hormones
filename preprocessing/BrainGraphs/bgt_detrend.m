function bgt_detrend(timeSeries, polyDegree, applyHPF, repetitionTime)
% Detrending of voxelwise timeseries.
% 
% FORMAT bgt_detrend(timeSeries, applyHPF, repetitionTime, polyDegree)
%
% REQUIRED INPUT:
%   timeSeries
%       String specifying file/path information for 4D (x,y,z + time)
%       image volume.
%
%   polyDegree
%       Degree of polynomial for detrending, usually = 1 (i.e. a linear
%       trend), but could also remove higher-order trends by editing the
%       value of this variable. Higher-order (e.g. quadratic) trends may
%       attenuate slow motion-related effects. However, this is generally
%       not recommended/necessary.
%
% OPTIONAL INPUT:
%   applyHPF
%       Determines whether or not to apply a 100s high-pass filter to the
%       timeseries. Enter 1 for YES.
%
%   repetitionTime
%       EPI sequence TR (entered in seconds). REQUIRED if applyHF = 1.
%
% OUTPUT
%   No explicit output. New images, along with detrending and high-pass
%   filtering matrices (if requested), will be written to the present
%   working directory.
%__________________________________________________________________________
%
% This function removes linear (or higher-order polynomial) trends from
% the timeseries. Useful for resting-state but not typically used when
% estimating task-related connectivity. It can also perform an initial
% high-pass filter using a discrete cosine series.
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

% Check input arguments.
%--------------------------------------------------------------------------

    if (nargin == 2), applyHPF = 0; repetitionTime = NaN; end
    if (nargin == 3), error('- ERROR: No repetitionTime specified! -'); end

% Load image volumes into memory, compress out non-brain voxels (implicit
% mask) for computational efficiency later.
%--------------------------------------------------------------------------

    V               = load_nii(timeSeries);
    data4D          = V.img;
    dataHeader      = V.hdr;
    [path,name,ext] = fileparts(timeSeries);
    
    [x,y,z,t]       = size(data4D);
    data2D          = reshape(data4D, x*y*z, t);
    mInd            = all(data2D == 0, 2);
    data2D(mInd,:)  = [];
    tsMatrix        = double(data2D)';
    
    clear V data2D data4D
    
    tsLength = size(tsMatrix,1);

    
% Apply initial high-pass filter via discrete cosine series if requested.
%--------------------------------------------------------------------------

    if applyHPF
        
        disp('|| Applying high-pass filter...');
            
        % Set filter cutoff to 100s (get rid of anything below 0.01 Hz).
            cutoff         = 100;
        
        % Total acquisition time.
            T              = tsLength*repetitionTime;
        
        % Determine minimum possible filter order.
            order          = floor((T/cutoff)*2)+1;
            
        % Demean the data.
            tsMu           = repmat(mean(tsMatrix), [size(txMatrix,1), 1]);
            tsMatrix       = tsMatrix - tsMu;
            
        % Create basis functions for filtering (pasted from spm_dctmtx, 
        % with thanks to Karl Friston).
            n              = (0:(tsLength-1))';
            C              = zeros(size(n,1),order);
            C(:,1)         = ones(size(n,1),1)/sqrt(tsLength);
            for k = 2:tsLength
                C(:,k) = sqrt(2/order)*cos(pi*(2*n+1)*(k-1)/(2*order));
            end

        % Filtering (i.e. 'residual-forming') matrix. Save this for later
        % just in case. 
            R              = eye(tsLength) - C*pinv(C);
            save('hpfMatrix','R');
            
        % Filtered timeseries with mean added back in.
            tsMatrix       = R'*tsMatrix + tsMu;
            
        disp('|| Finished filtering');    
        
        clear repetitionTime cuttoff T order n k C R
            
    end

    
% Remove polynomial (typically linear) trend from voxelwise timeseries.
%--------------------------------------------------------------------------

    disp('|| Detrending voxelwise timeseries...');
        
    % Linear basis function(s).
        linearBasis = repmat((1:tsLength)', [1, polyDegree]);
        
    % Polynomial contrast weights.
        polyWeight  = repmat(1:polyDegree, [tsLength 1]);
        
    % Final weighted regressors: m x n matrix where m = tsLength and
    % n = 1 + polyDegree. First column will always be constant (to remove 
    % mean), remaining column(s) will contain linear/polynomial basis.
        C           = [ones(tsLength,1), linearBasis.^polyWeight];
        
    % Detrending matrix - save for later.
        R           = eye(tsLength) - C*pinv(C);
        save('dtMatrix','R');
            
    % Detrended timeseries.
        tsMatrix    = R'*tsMatrix;
        
    disp('|| Finished detrending');
        
    clear polyDegree linearBasis polyWeight C R
    
    
% Write output.
%--------------------------------------------------------------------------

    data2D              = zeros(t, x*y*z);
    data2D(:,mInd == 0) = tsMatrix;
    data4D              = reshape(data2D', [x y z t]);
            
    nii.img             = data4D;
    nii.hdr             = dataHeader;
                
    save_nii(nii, [path '/d' name ext]);
                
end