function bgt_bandpass(timeSeries, repetitionTime, lowerFreq, upperFreq)
% Temporal bandpass filtering of a timeseries.
%
% FORMAT bgt_bandpass(timeSeries, repetitionTime, lowerFreq, upperFreq)
%
% REQUIRED INPUT:
%   timeSeries
%       Data can be specified here in one of two forms:
%           1) String specifying file/path information for 4D
%              (x,y,z + time) image volume.
%              NOTE: must be a character array, not a cell array.
%           2) m x n matrix where m = number of TRs and n = number of
%              regions of interest (ROIs) in your network.
%
%   repetitionTime
%       EPI sequence TR (entered in seconds).
%
%   lowerFreq
%       Highpass filter cutoff (entered in Hz - a common value is .01).
%
%   upperFreq
%       Lowpass filter cutoff (entered in Hz - a common value is .10).
%
% OUTPUT:
%   No explicit output. If timeSeries was entered as an image volume, new
%   images will be written to the working directory; if timeSeries was
%   entered as a matrix, a new .mat file will be written to the working
%   directory.
%__________________________________________________________________________
%
% This function will apply a temporal bandpass filter to the timeseries.
% This serves as an alternative to the maximal overlap discrete wavelet
% transform (MODWT), which may offer extra denoising but is limited to
% certain frequency octaves (depending on the TR). Here, any range of
% frequencies can be extracted as specified in the function arguments.
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
                t        = size(tsMatrix, 1);
                
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
    
    
% Apply temporal filter given the frequency parameters defined above.
%--------------------------------------------------------------------------

    disp(['|| Applying temporal bandpass filter: ' num2str(lowerFreq) ' - ' num2str(upperFreq) ' Hz...']);

    % First demean the data.
        tsMu     = repmat(mean(tsMatrix), [size(tsMatrix,1), 1]);
        tsMatrix = tsMatrix - tsMu;
        
    % Zero-pad for faster FFT.
        targetLength = 2^nextpow2(t);
        tsMatrix     = [tsMatrix; zeros(targetLength - t, size(tsMatrix,2))];
        
    % Get FFT, apply filter model, then invert.
        tsFFT           = fft(tsMatrix, [], 1);
        tIndex          = 0:size(tsFFT,1)-1;
        f               = min(tIndex, size(tsFFT,1)-tIndex);
        fIndex          = find(f < lowerFreq*(repetitionTime*size(tsFFT,1)) | f > upperFreq*(repetitionTime*size(tsFFT,1)));
        fIndex          = fIndex(fIndex > 1);
        tsFFT(fIndex,:) = 0;
        bandpassSeries  = real(ifft(tsFFT, [], 1));
        
    % Curtail end of zero-padded timeseries, add back mean.
        bandpassSeries = bandpassSeries(1:t,:) + tsMu;        
        
    disp('|| Finished filtering');    
        
    
% Write output.
%--------------------------------------------------------------------------

    switch tsType
        
        case 0
            
            save('bandpassSeries','bandpassSeries');
            
        case 1
                        
            data2D              = zeros(t, x*y*z);
            data2D(:,mInd == 0) = bandpassSeries;
            data4D              = reshape(data2D', [x y z t]);
            
            nii.img             = data4D;
            nii.hdr             = dataHeader;
            
            save_nii(nii, [path '/b' name ext]);
            
    end       

end