function jacobs28andMe_processData(whichDays)
% Run through standard preprocessing routine.
%
% FORMAT jacobs28andMe_processData(whichDays)
%
% REQUIRED INPUT:
%   whichDays
%       Cell array of strings defining which sessions to preprocess.
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
        
        % Apply the following preprocessing steps for this session:
        
        jacobs28andMe_calculateVDM;
        jacobs28andMe_realignUnwarp;
        jacobs28andMe_segmentStrip;
        jacobs28andMe_estimateAnatomicalNoise;
        jacobs28andMe_coregisterTemplate;
        jacobs28andMe_smooth;
        
        % Jump back to the parent directory.
            
        cd(parentDir);
                    
    end
    
% Display total computation time.
%--------------------------------------------------------------------------

    procEnd = toc(procStart);
    disp(['|| Jobs completed for ' num2str(nSessions) ' sessions in ' num2str(procEnd/60) ' minutes']);
    
end