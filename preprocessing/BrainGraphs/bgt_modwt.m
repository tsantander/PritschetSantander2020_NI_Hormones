function bgt_modwt(timeSeries, repetitionTime, wavScale)
% Temporal filtering of a timeseries via MODWT.
%
% FORMAT bgt_modwt(timeSeries, repetitionTime, wavScale)
%
% REQUIRED INPUT:
%   timeSeries
%       Data can be specified here in one of two forms:
%           1) String specifying file/path information for 4D
%              (x,y,z + time) image volume.
%              NOTE: must be a character array, not a cell array.
%           2) m x n matrix where m = number of TRs and n = number of
%              nodes (i.e. ROIs) in your network.
%
%   repetitionTime
%       EPI sequence TR (entered in seconds).
%
%   wavScale
%       Frequency band (i.e. wavelet scale) to extract. This argument can
%       be entered in one of two ways:
%           1) Scalar value (e.g. 2), which will decompose the timeseries
%              down to the octave specified. The resulting data will be a
%              matrix (or 4D image) of WAVELET COEFFICIENTS.
%           2) Vector (e.g. 2:3), which will decompose the timeseries to
%              the maximum octave specified - this will then perform an
%              inverse transform to RECONSTRUCT THE ORIGINAL SIGNAL with
%              the unspecified wavlet scales removed (e.g. scale 1 if 2:3
%              was given). In theory, this could allow you to estimate
%              connectivity over a wider frequency band. While it is most 
%              conventional to use wavelet coefficients from a single
%              scale, this can be highly restrictive when data were
%              acquired at high temporal resolution (e.g. with multiband
%              EPI).
%
% OUTPUT:
%   No explicit output. If timeSeries was entered as an image volume, new
%   images will be written to the working directory; if timeSeries was
%   entered as a matrix, a new .mat file will be written to the working
%   directory.
%__________________________________________________________________________
%
% This function will filter a timeseries using a maximal overlap discrete
% wavelet transform (MODWT). Signal decomposition via MODWT can be
% particularly well suited for nonstationary timeseries such as fMRI
% because it offers increased sensitivity to small signal changes in a
% noisy background. However, note that the frequency bands we can extract
% are limited by our sampling rate - you can use the helper function
% bgt_wavCalc to determine which octaves are obtainable given the TR of
% your pulse sequence.
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
                
                error('- ERROR: timeSeries argument must be a character array or matrix -');
                
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
            
            error('- ERROR: timeSeries argument must be a character array or matrix -');
    
    end
    
    clear timeSeries
    
% Perform wavelet decomposition for the scale(s) specified - return wavelet
% coefficients for a single scale, or reconstruct the signal within a range
% if necessary.
%--------------------------------------------------------------------------

    allScales = zeros(max(wavScale),2);

    for iWav = 1:max(wavScale)
            
        allScales(iWav,:) = [(1/(2^(iWav+1)*repetitionTime)), (1/(2^(iWav)*repetitionTime))];
        
    end 
    
    usedRange = allScales(wavScale,:); clear allScales
    lower     = min(reshape(usedRange,1,numel(usedRange)));
    upper     = max(reshape(usedRange,1,numel(usedRange)));
    
    % Compute wavelet decompositions using MODWT: note that we use a 
    % Daubechies filter with length = 8, per Zhang et al (2016, PLoS One).
    
        disp(['|| Extracting frequencies: ~' num2str(round(lower,3)) '-' num2str(round(upper,3)) ' Hz...']);
    
        [WJt, VJt, att] = modwt(tsMatrix, 'd8', max(wavScale), 'circular', 'RetainVJ', 1);
    
    % Get relevant data based on wavScale.
    
        if (isscalar(wavScale))
            
            waveletSeries = squeeze(WJt(:,wavScale,:));
            
        else
            
            wtfS = modwt_filter(att.WTF);
            ht   = wtfS.h;
            gt   = wtfS.g;
            
            waveletSeries = zeros(size(tsMatrix));
            
            for iChannel = 1:size(tsMatrix,2)
                
                Vin = VJt(:,max(wavScale),iChannel);
                
                % Reconstruct starting with lowest-frequency harmonics.
                
                for jScale = max(wavScale):-1:min(wavScale)
                    
                    Vout = imodwtj(WJt(:,jScale,iChannel), Vin, ht, gt, jScale);
                    Vin  = Vout;
                    
                end
                
                waveletSeries(:,iChannel) = Vout;
        
            end
            
        end
        
    clear WJt VJt att
    
    disp('|| Wavelet decomposition complete');
        
% Write output.
%--------------------------------------------------------------------------

    switch tsType
        
        case 0
            
            save('waveletSeries','waveletSeries');
            
        case 1
                        
            data2D              = zeros(t, x*y*z);
            data2D(:,mInd == 0) = waveletSeries;
            data4D              = reshape(data2D', [x y z t]);
            
            nii.img             = data4D;
            nii.hdr             = dataHeader;
            
            save_nii(nii, [path '/t' name ext]);
            
            
    end        

end