function jacobs28andMe_restingState(whichDays)
% Run additional preprocessing steps for rs-fMRI and estimate coherence.
%
% FORMAT jacobs28andMe_restingState(whichDays)
%
% REQUIRED INPUT:
%   whichDays
%       Cell array of strings defining which sessions to run.
%__________________________________________________________________________
%
% Author:
%   Tyler Santander (t.santander@psych.ucsb.edu)
%   Institute for Collaborative Biotechnologies
%   Department of Psychological & Brain Sciences
%   University of California, Santa Barbara
%   August 2018
%__________________________________________________________________________

% How many days worth of data are we trying to cover?
%--------------------------------------------------------------------------

    nSessions = length(whichDays);
    
% Start the clock so we can track overall computation time.
%--------------------------------------------------------------------------

    procStart = tic;
    
% Begin looping over sessions...
%--------------------------------------------------------------------------
    
    disp(' ');
    
    parentDir = pwd;
    
    for iDay = 1:nSessions
        
        dayID = whichDays{iDay};
        
        disp(['|| Running day: ' dayID '. Please wait...']);
        
        cd(dayID);
        
        restingDir = [pwd '/data.functional.rest'];
        
        % Create a new directory in which to store the results, copy over
        % data and navigate into it.
        
        mkdir('results.network.rest');
        cd('results.network.rest');
        
        unix(['mv ' restingDir '/' spm_select('List', restingDir, '^sruf.*nii.gz$') ' ' pwd '/rest4D.nii.gz']);
        unix(['cp ' restingDir '/rp* ' pwd]);
        unix(['cp ' restingDir '/anat* ' pwd]);
        
        % Scale functional data to grand median of 1000.
        
        bgt_globalNorm([pwd '/' spm_select('List', pwd, '^rest4D.nii.gz$')], 'median');
        
        % Linearly detrend voxelwise timeseries.
        
        bgt_detrend([pwd '/' spm_select('List', pwd, '^grest4D.nii.gz$')], 1);
        
        % Detrend motion parameters / anatomical noise and regress from
        % timeseries.
        
        load('dtMatrix.mat' ,'R');
        load('anatNoise.mat', 'anatNoise');
        
        motionParams = load([pwd '/' spm_select('List', pwd, '^rp.*txt$')]);
        dtMotion     = R'*motionParams;
        dtAnatNoise  = R'*anatNoise;
        
        bgt_regressNuisance([pwd '/' spm_select('List', pwd, '^dgrest4D.nii.gz$')], dtMotion, 'fristonAR1', dtAnatNoise);
        
        % Extract regional timeseries.
        
        [timeSeries] = bgt_extractRegionalTimeseries([pwd '/' spm_select('List', pwd, '^ndgrest4D.nii.gz$')], ...
                        ['/home/tyler/28andMe/lpTemplate/atlases/' spm_select('List', '/home/tyler/28andMe/lpTemplate/atlases', '^compositeAtlasWarp.*nii.gz$')], ...
                        'eigen1');
                    
        % Get relevant frequency band for modwt and decompose.
        
        wavScales  = 3:6;
        [wavFreqs] = bgt_wavCalc(.720, 6);
        freqBand   = [wavFreqs(6,2), wavFreqs(3,3)];
                
        bgt_modwt(timeSeries, .720, wavScales);
        
        load([pwd '/waveletSeries.mat'], 'waveletSeries');
        
        % Estimate coherence and apply FDR correction.
        
        [coherence] = bgt_coherenceMatrix(waveletSeries, .720, freqBand, 'Welch', 'parametric', 0);
        bgt_coherenceFDR(coherence, .05);
        
        % Navigate back to parent directory.
                            
        cd(parentDir); disp(' ');
                    
    end
    
% Display total computation time.
%--------------------------------------------------------------------------

    procEnd = toc(procStart);
    disp(['|| Jobs completed for ' num2str(nSessions) ' sessions in ' num2str(procEnd/60) ' minutes']);
    
end