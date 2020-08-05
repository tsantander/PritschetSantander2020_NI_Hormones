function jacobs28andMe_smooth
% Apply spatial smoothing.
%
% FORMAT jacobs28andMe_smooth
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
    
% Get all the realigned/unwarped functional data in 2mm template space.
%--------------------------------------------------------------------------

    restEPI = cellstr(strcat([pwd '/data.functional.rest/'], spm_select('List', [pwd '/data.functional.rest'], '^rufrest4D.*nii$')));

% Specify smoothing kernel (4mm FWHM) and run.
%--------------------------------------------------------------------------

    matlabbatch{1}.spm.spatial.smooth.data   = restEPI;
    matlabbatch{1}.spm.spatial.smooth.fwhm   = [4 4 4];
    matlabbatch{1}.spm.spatial.smooth.dtype  = 0;
    matlabbatch{1}.spm.spatial.smooth.im     = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = 's';
    
    spm_jobman('run',matlabbatch);
    
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
