function [coherence] = bgt_coherenceMatrix(filteredSeries, repetitionTime, freqBand, cohMethod, sigTest, fig)
% Compute functional association via magnitude-squared coherence.
%
% FORMAT [coherence] = bgt_coherenceMatrix(filteredSeries, repetitionTime, freqBand, cohMethod, sigTest, fig)
%
% REQUIRED INPUT:
%   filteredSeries
%       m x n matrix of temporally-filtered data, where m = number of
%       timepoints and n = number of regions of interest (ROIs) in the 
%       network.
%
%   repetitionTime
%       EPI sequence TR (entered in seconds).
%
%   freqBand
%       Frequency band over which to estimate coherence. This should be a
%       vector (e.g. [0.01, 0.10]).
%
%   cohMethod
%       Method used to compute functional coherence between ROIs. This is
%       ultimately a matter of personal preference. Specify (in single
%       quotes) either:
%           'Welch'             -   Matlab's mscohere function using
%                                   Welch's averaged modified periodogram
%                                   method.
%           'MVDR'              -   Minimum variance distortionless
%                                   response (Benesty, Chen, & Huang,
%                                   2006). Like Welch's, this is a
%                                   nonparametric estimation algorithm;
%                                   however, the bandpass filters of the
%                                   MVDR spectrum are both data-  and
%                                   frequency-dependent, in contrast to
%                                   Welch's method, which instead employs
%                                   data- and frequency-independent
%                                   discrete Fourier matrices. This
%                                   technique may offer higher resolution
%                                   and more reliable results, but it is
%                                   not yet commonly used in coherence
%                                   estimation for fMRI. NB: if you choose
%                                   MVDR, you MUST perform significance
%                                   testing via surrogates.
%
%   sigTest
%       Indicates a method for calculating p-values (i.e. the probability
%       of obtaining a given coherence value, assuming the null of no
%       association between brain regions). Broadly, these fall into
%       parametric or surrogate (nonparametric) approaches. 'Parametric'
%       estimation is fastest. 'Surrogate' estimation (i.e. nonparametric 
%       permutation testing or amplitude-adjusted Fourier transforms) 
%       allows us to obtain empirical null distributions, which I generally
%       prefer, but it comes at computational cost. NOTE: surrogate testing
%       is REQUIRED for coherence estimation via MVDR! Specify (in single 
%       quotes):
%           'parametric'        -   Uses formulas taken from Shumway &
%                                   Stoffer (2010) to estimate the
%                                   theoretical p-values associated with
%                                   a wide range of coherence values.
%                                   Based on the F-distribution.
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
%                                   linear structure of the signal
%                                   (important because coherence is a
%                                   linear estimate of association).
%           'IAAFT'             -   Uses iterative AAFT to generate
%                                   surrogate data. This will take longer
%                                   than standard AAFT but provides a
%                                   'better' approximation of the original
%                                   signal's autocorrelation function in
%                                   addition to its amplitude distribution.
%
%   fig
%       Indicates whether or not to display the coherence heatmap. Enter
%       1 for YES or 0 for NO.
%
% OUTPUT:
%   coherence
%       Structure array with the following fields:
%           .uncorrected
%               .associationMatrix     -    Symmetric n x n matrix whose
%                                           elements (i,j) indicate the
%                                           coherence between the ith and
%                                           jth ROIs in your network.
%               .pValues               -    Symmetric n x n matrix whose
%                                           elements (i,j) indicate the
%                                           probability of obtaining the
%                                           coherence value contained in
%                                           associationMatrix(i,j) by
%                                           chance (if the null of no
%                                           association is true).
%__________________________________________________________________________
%
% This function will compute the extent of co-activity between pairs of
% brain regions using magnitude-squared coherence. This approach may be
% preferable to standard Pearson product-moment correlations because it:
%   1) Allows us to compute frequency-specific covariances.
%   2) Is normalized to [0,1].
%   3) Is robust to differences in the hrf between brain regions.
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

% Initial setup.
%--------------------------------------------------------------------------

    nROI  = size(filteredSeries,2);
    nEdge = (nROI^2 - nROI)/2;
    lower = freqBand(1);
    upper = freqBand(2);

    
% Specify parameters for coherence computations - this will depend on the
% algorithm chosen: Matlab's default (Welch) vs. MVDR.
%--------------------------------------------------------------------------

    % FFT length - determines the freqs at which coherence is estimated.
    % This is required regardless of the algorithm specified. There are 
    % data-driven means of identifying the optimal nfft, but a standard
    % default is one that best approximates 8 equal sections of our signal
    % to the nearest power of 2. In order to ensure we're able to do this
    % we'll first zero-pad the timeseries to the nearest power of 2 (this
    % will also help speed up FFT computations later).
        
        tsLength       = size(filteredSeries,1);
        
        targetLength   = 2^nextpow2(tsLength);
        filteredSeries = [filteredSeries; zeros(targetLength - tsLength, nROI)];
                        
        nfft           = 2^round(log2(targetLength/8));
        
    % Now specify additional parameters based on the algorithm of choice.
    
        switch cohMethod
            
            case 'Welch'
                
                % Compute a periodic Hann window.
            
                    nfftPlus    = nfft + 1;
                    approxHalf  = ceil(nfftPlus/2);
                    window      = (0:approxHalf-1)'/(nfftPlus-1);
                    window      = 0.5 - 0.5*cos(2*pi*window);
                    window      = [window; window(end-1:-1:1)];
                    window(end) = [];
                    
                % Number of samples by which the sections overlap (fairly
                % standard to use 50%, so nfft/2).
    
                    nOverlap   = nfft/2;
                    
                % Sampling frequency.
            
                    sampleFreq = 1/repetitionTime;
                    
            case 'MVDR'
                
                % Resolution parameter (in theory, the higher the better,
                % but we'll just take one order of magnitude larger than
                % our current nfft). We'll also compute the corresponding
                % frequencies at this resolution while we're at it.
                
                    resolution = 10*nfft;
                    [freq]     = 0:1/resolution:1-1/resolution;
                    
        end
        
    % If using a parametric test of significance, determine cutoff scores
    % for a range of p-values.
    
        if ((strcmp(cohMethod, 'Welch') == 1) && (strcmp(sigTest, 'parametric') == 1))
                            
            alpha   = .0001:.0001:.99;                
            df      = 2*nOverlap*tsLength/targetLength;
            critF   = finv(1-alpha, 2, df-2);
            critC   = critF./(df/2-1+critF);
            
        elseif ((strcmp(cohMethod, 'Welch') == 0) && (strcmp(sigTest, 'parametric') == 1))
            
            error('- ERROR: Parametric significance testing not applicable for MVDR! -');
                                        
        end
        
    % Initialize arrays to store results.

        coherence         = [];
        associationMatrix = zeros(nROI);
        pValues           = zeros(nROI);
        
        % Estimate coherence between nodal timeseries.
        
            startCoh    = tic;
            edgeCounter = 1;
            
            disp(['|| Computing magnitude-squared coherence at ~' num2str(round(lower,3)) '-' num2str(round(upper,3)) ' Hz...']);
            
            textprogressbar(['|| Estimated progress over ' num2str(nEdge) ' network edges: ']);

            for iROI = 1:nROI
                
                for jROI = 1:nROI
                    
                    if (iROI <= jROI)   % Skip diagonal and lower triangle.
                        
                        continue
                        
                    else                % Otherwise, estimate coherence.
                        
                        progress = 100*(edgeCounter/nEdge);
                        textprogressbar(progress);
                    
                        switch cohMethod
                        
                            case 'Welch'
                            
                                [coh, freq] = mscohere(filteredSeries(:,iROI), filteredSeries(:,jROI), window, nOverlap, nfft, sampleFreq);
                            
                            case 'MVDR'
                            
                                [coh] = coherence_MVDR(filteredSeries(:,iROI), filteredSeries(:,jROI), nfft, resolution);
                            
                        end
                        
                        % Get average coherence within our band of
                        % interest.
                    
                        lowerID      = find(freq >= lower, 1, 'first');
                        upperID      = find(freq <= upper, 1, 'last');
                        avgCoherence = mean(coh(lowerID:upperID));
                    
                        associationMatrix(iROI,jROI) = avgCoherence;
                        
                        % Get a p-value.
                        
                        switch sigTest
                            
                            case 'parametric'
                    
                                try
                                    pValues(iROI,jROI) = alpha(find(critC < avgCoherence, 1, 'first'));
                                catch
                                    pValues(iROI,jROI) = 1;
                                end
                                
                            case {'NPT', 'AAFT', 'IAAFT'}
                                
                                % Generate data to empirically test null
                                % distributions of coherence. If 'NPT',
                                % we'll use the 'jROI' timeseries and
                                % randomly permute it (doesn't really
                                % matter because coherence is symmetric);
                                % if 'AAFT' or 'IAAFT' we need to use both
                                % because we're generating surrogate data.
                                
                                nIter   = 1000;
                                randCoh = zeros(nIter,1);
                                
                                switch sigTest
                                    
                                    case 'NPT'
                                        
                                        temp  = filteredSeries(1:tsLength,jROI);
                                        yPerm = zeros(tsLength, nIter);
                                        
                                        for iPerm = 1:nIter
                                            yPerm(:,iPerm) = temp(randperm(tsLength));
                                        end
                                        
                                        yPerm = [yPerm; zeros(targetLength - tsLength, nIter)];
                                        
                                        clear temp
                                        
                                        % Loop over permutations.
                                
                                        for iPerm = 1:nIter
                                                                    
                                            switch cohMethod
                        
                                                case 'Welch'
                            
                                                    [coh, freq] = mscohere(filteredSeries(:,iROI), yPerm(:,iPerm), window, nOverlap, nfft, sampleFreq);
                            
                                                case 'MVDR'
                            
                                                    [coh] = coherence_MVDR(filteredSeries(:,iROI), yPerm(:,iPerm), nfft, resolution);
                                            
                                            end
                                    
                                            lowerID        = find(freq >= lower, 1, 'first');
                                            upperID        = find(freq <= upper, 1, 'last');                    
                                            randCoh(iPerm) = mean(coh(lowerID:upperID));
                                    
                                        end
                                        
                                        clear yPerm
                                            
                                    case 'AAFT'
                                        
                                        xPerm = AAFT(filteredSeries(1:tsLength,iROI), nIter);
                                        yPerm = AAFT(filteredSeries(1:tsLength,jROI), nIter);
                                        
                                        xPerm = [xPerm; zeros(targetLength - tsLength, nIter)];
                                        yPerm = [yPerm; zeros(targetLength - tsLength, nIter)];
                                        
                                        % Loop over surrogates.
                                
                                        for iPerm = 1:nIter
                                                                    
                                            switch cohMethod
                        
                                                case 'Welch'
                            
                                                    [coh, freq] = mscohere(xPerm(:,iPerm), yPerm(:,iPerm), window, nOverlap, nfft, sampleFreq);
                            
                                                case 'MVDR'
                            
                                                    [coh] = coherence_MVDR(xPerm(:,iPerm), yPerm(:,iPerm), nfft, resolution);
                                            
                                            end
                                    
                                            lowerID        = find(freq >= lower, 1, 'first');
                                            upperID        = find(freq <= upper, 1, 'last');                    
                                            randCoh(iPerm) = mean(coh(lowerID:upperID));
                                            
                                        end
                                        
                                        clear xPerm yPerm
                                        
                                    case 'IAAFT'
                                        
                                        xPerm = IAAFT(filteredSeries(1:tsLength,iROI), nIter);
                                        yPerm = IAAFT(filteredSeries(1:tsLength,jROI), nIter);
                                        
                                        xPerm = [xPerm; zeros(targetLength - tsLength, nIter)];
                                        yPerm = [yPerm; zeros(targetLength - tsLength, nIter)];
                                        
                                        % Loop over surrogates.
                                
                                        for iPerm = 1:nIter
                                                                    
                                            switch cohMethod
                        
                                                case 'Welch'
                            
                                                    [coh, freq] = mscohere(xPerm(:,iPerm), yPerm(:,iPerm), window, nOverlap, nfft, sampleFreq);
                            
                                                case 'MVDR'
                            
                                                    [coh] = coherence_MVDR(xPerm(:,iPerm), yPerm(:,iPerm), nfft, resolution);
                                            
                                            end
                                    
                                            lowerID        = find(freq >= lower, 1, 'first');
                                            upperID        = find(freq <= upper, 1, 'last');                    
                                            randCoh(iPerm) = mean(coh(lowerID:upperID));
                                            
                                        end
                                        
                                        clear xPerm yPerm
                                        
                                end
                                
                                % Get an empirical p-value based on the
                                % proportion of permuted coherences that 
                                % were greater than or equal to the 'true' 
                                % estimate for this pair of timeseries.
                                
                                pValues(iROI,jROI) = sum(double(randCoh >= avgCoherence))/nIter;
                            
                        end
                        
                        edgeCounter = edgeCounter + 1;
                        
                    end
                
                end
    
            end
            
            textprogressbar(' DONE');
            
            endCoh = toc(startCoh);
            disp(['|| Coherence between ' num2str(nEdge) ' network edges computed in ' num2str(round(endCoh/60,2)) ' minutes']);
            
            coherence.uncorrected.associationMatrix = associationMatrix + associationMatrix';
            coherence.uncorrected.pValues           = pValues + pValues';
                        
    % Display coherence heatmap (optional, set 'fig' argument to 1 if YES).
    
        if fig
            
            figure; image(coherence.uncorrected.associationMatrix .* 64);
            
        end

        
% Save coherence structure.
%--------------------------------------------------------------------------

    save('coherence', 'coherence');

end