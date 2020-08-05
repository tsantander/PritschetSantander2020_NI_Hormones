function jacobs28andMe_realignUnwarp
% Perform motion realignment and unwarp EPIs using the VDM.
%
% FORMAT jacobs28andMe_realignUnwarp
%__________________________________________________________________________
%
% Author:
%   Tyler Santander (t.santander@psych.ucsb.edu)
%   Institute for Collaborative Biotechnologies
%   Department of Psychological & Brain Sciences
%   University of California, Santa Barbara
%   August 2018
%__________________________________________________________________________

% Initialize default SPM configurations for fMRI.
%--------------------------------------------------------------------------

    setDefaultsSPM;
    
    matlabbatch = {};
        
% Tell SPM where the scans live (in order of acquisition).
%--------------------------------------------------------------------------

    % First let's grab the voxel displacement map since we're using the 
    % same one for all of our EPIs.
    
        vdm = cellstr(strcat([pwd '/data.fieldmap2/'], ...
                spm_select('List', [pwd '/data.fieldmap2'], '^vdm.*nii$')));

    % Now resting-state.
    
        matlabbatch{1}.spm.spatial.realignunwarp.data(1).scans  = cellstr(strcat([pwd '/data.functional.rest/'], ...
                                                                    spm_select('List', [pwd '/data.functional.rest'], '^f.*nii$')));
        matlabbatch{1}.spm.spatial.realignunwarp.data(1).pmscan = vdm;
        
    % And the remaining EPIs.
    
        matlabbatch{1}.spm.spatial.realignunwarp.data(2).scans  = cellstr(strcat([pwd '/data.functional.attention1/'], ...
                                                                    spm_select('List', [pwd '/data.functional.attention1'], '^f.*nii$')));
        matlabbatch{1}.spm.spatial.realignunwarp.data(2).pmscan = vdm;
        
        matlabbatch{1}.spm.spatial.realignunwarp.data(3).scans  = cellstr(strcat([pwd '/data.functional.attention2/'], ...
                                                                    spm_select('List', [pwd '/data.functional.attention2'], '^f.*nii$')));
        matlabbatch{1}.spm.spatial.realignunwarp.data(3).pmscan = vdm;
        
        matlabbatch{1}.spm.spatial.realignunwarp.data(4).scans  = cellstr(strcat([pwd '/data.functional.attention3/'], ...
                                                                    spm_select('List', [pwd '/data.functional.attention3'], '^f.*nii$')));
        matlabbatch{1}.spm.spatial.realignunwarp.data(4).pmscan = vdm;
        
        matlabbatch{1}.spm.spatial.realignunwarp.data(5).scans  = cellstr(strcat([pwd '/data.functional.attention4/'], ...
                                                                    spm_select('List', [pwd '/data.functional.attention4'], '^f.*nii$')));
        matlabbatch{1}.spm.spatial.realignunwarp.data(5).pmscan = vdm;
        
        matlabbatch{1}.spm.spatial.realignunwarp.data(6).scans  = cellstr(strcat([pwd '/data.functional.reward/'], ...
                                                                    spm_select('List', [pwd '/data.functional.reward'], '^f.*nii$')));
        matlabbatch{1}.spm.spatial.realignunwarp.data(6).pmscan = vdm;
        
        if (exist([pwd '/data.functional.criterion1'], 'dir') == 7)
        
            matlabbatch{1}.spm.spatial.realignunwarp.data(7).scans  = cellstr(strcat([pwd '/data.functional.criterion1/'], ...
                                                                        spm_select('List', [pwd '/data.functional.criterion1'], '^f.*nii$')));
            matlabbatch{1}.spm.spatial.realignunwarp.data(7).pmscan = vdm;
        
            matlabbatch{1}.spm.spatial.realignunwarp.data(8).scans  = cellstr(strcat([pwd '/data.functional.criterion2/'], ...
                                                                        spm_select('List', [pwd '/data.functional.criterion2'], '^f.*nii$')));
            matlabbatch{1}.spm.spatial.realignunwarp.data(8).pmscan = vdm;
            
        end
    
% Define all parameters for realign/unwarp estimation.
%--------------------------------------------------------------------------

    % Motion estimation for realignment.

        matlabbatch{1}.spm.spatial.realignunwarp.eoptions.quality = 1;
        matlabbatch{1}.spm.spatial.realignunwarp.eoptions.sep     = 4;
        matlabbatch{1}.spm.spatial.realignunwarp.eoptions.fwhm    = 5;
        matlabbatch{1}.spm.spatial.realignunwarp.eoptions.rtm     = 0;
        matlabbatch{1}.spm.spatial.realignunwarp.eoptions.einterp = 7;
        matlabbatch{1}.spm.spatial.realignunwarp.eoptions.ewrap   = [0 0 0];
        matlabbatch{1}.spm.spatial.realignunwarp.eoptions.weight  = '';

    % Deformation field estimation for unwarp.

        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.basfcn   = [12 12];
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.regorder = 1;
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.lambda   = 100000;
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.jm       = 0;
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.fot      = [4 5];
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.sot      = [];
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.uwfwhm   = 4;
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.rem      = 1;
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.noi      = 5;
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.expround = 'Average';

    % Interpolation/reslicing.

        matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.uwwhich = [2 1];
        matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.rinterp = 7;
        matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.wrap    = [0 0 0];
        matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.mask    = 1;
        matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.prefix  = 'u';

% Run motion/distortion correction.
%--------------------------------------------------------------------------

    spm_jobman('run', matlabbatch);

% Create new folder for the mean functional image and move it there.
%--------------------------------------------------------------------------

    mkdir('./data.functional.mean');
    unix('mv ./data.functional.rest/meanuf* ./data.functional.mean');

% Estimate framewise displacement for each scanning run.
%--------------------------------------------------------------------------

    % Initialize structure for resting-state.

        framewiseDisplacement.rest = [];
        
    % Obtain realignment parameters.    
    
        rpData = load(strcat([pwd '/data.functional.rest/'], ...
                        spm_select('List', [pwd '/data.functional.rest'], '^rp.*txt$')));
        
    % Compute and store. 
        
        rpData(:,4:6) = rpData(:,4:6) .* 50;
        dx            = diff(rpData);
        fwd           = sum(abs(dx),2);
            
        framewiseDisplacement.rest.series = fwd;
        framewiseDisplacement.rest.mean   = mean(fwd);
        framewiseDisplacement.rest.max    = max(fwd);

    % Loop over the selective attention runs.

        for iRun = 1:4
        
            % Initialize structure for this run.

                framewiseDisplacement.(['attention' num2str(iRun)]) = [];
        
            % Obtain realignment parameters.    
    
                rpData = load(strcat([pwd '/data.functional.attention' num2str(iRun) '/'], ...
                                spm_select('List', [pwd '/data.functional.attention' num2str(iRun)], '^rp.*txt$')));
        
            % Compute and store. 
        
                rpData(:,4:6) = rpData(:,4:6) .* 50;
                dx            = diff(rpData);
                fwd           = sum(abs(dx),2);
            
                framewiseDisplacement.(['attention' num2str(iRun)]).series = fwd;
                framewiseDisplacement.(['attention' num2str(iRun)]).mean   = mean(fwd);
                framewiseDisplacement.(['attention' num2str(iRun)]).max    = max(fwd);
            
        end
        
    % Reward localizer.
    
        framewiseDisplacement.reward = [];
            
        rpData = load(strcat([pwd '/data.functional.reward/'], ...
                        spm_select('List', [pwd '/data.functional.reward'], '^rp.*txt$')));
                
        rpData(:,4:6) = rpData(:,4:6) .* 50;
        dx            = diff(rpData);
        fwd           = sum(abs(dx),2);
            
        framewiseDisplacement.reward.series = fwd;
        framewiseDisplacement.reward.mean   = mean(fwd);
        framewiseDisplacement.reward.max    = max(fwd);
        
    % Criterion runs.
    
        if (exist([pwd '/data.functional.criterion1'], 'dir') == 7)

            for iRun = 1:2
        
                % Initialize structure for this run.

                    framewiseDisplacement.(['criterion' num2str(iRun)]) = [];
        
                % Obtain realignment parameters.    
    
                    rpData = load(strcat([pwd '/data.functional.criterion' num2str(iRun) '/'], ...
                                    spm_select('List', [pwd '/data.functional.criterion' num2str(iRun)], '^rp.*txt$')));
        
                % Compute and store. 
        
                    rpData(:,4:6) = rpData(:,4:6) .* 50;
                    dx            = diff(rpData);
                    fwd           = sum(abs(dx),2);
            
                    framewiseDisplacement.(['criterion' num2str(iRun)]).series = fwd;
                    framewiseDisplacement.(['criterion' num2str(iRun)]).mean   = mean(fwd);
                    framewiseDisplacement.(['criterion' num2str(iRun)]).max    = max(fwd);
                    
            end
            
        end
            
    % Save structure.
        
        mkdir('./data.fwd');
        cd('./data.fwd');
        save framewiseDisplacement framewiseDisplacement;
        cd ..

end

%-------------------------------------------------------------------------%
% BEGIN SUBROUTINES                                                       %
%-------------------------------------------------------------------------%

% Initialize default parameters for SPM.
%-------------------------------------------------------------------------%
function setDefaultsSPM

    spm('defaults','fMRI');
    warning off MATLAB:FINITE:obsoleteFunction;
    spm_jobman('initcfg');
    
end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% END SUBROUTINES                                                         %
%-------------------------------------------------------------------------%