function jacobs28andMe_coregisterTemplate
% Register T1 and EPIs to subject-specific anatomical template space.
%
% FORMAT jacobs28andMe_coregisterTemplate
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
    
% Specify N3-regularised mean EPI as the REFERENCE SCAN.
%--------------------------------------------------------------------------
    
    matlabbatch{1}.spm.spatial.coreg.estimate.ref = cellstr(strcat([pwd '/data.functional.mean/'], spm_select('List', [pwd '/data.functional.mean'], '^meanN3.*nii$')));
    
% Specify skull-stripped hires as SOURCE SCAN.
%--------------------------------------------------------------------------

    matlabbatch{1}.spm.spatial.coreg.estimate.source = cellstr(strcat([pwd '/data.anatomical.hires/'], spm_select('List', [pwd '/data.anatomical.hires'], '^ss-ms.*nii$')));
    
    matlabbatch{1}.spm.spatial.coreg.estimate.other             = {''};
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep      = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol      = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm     = [7 7];
    
% Write transformation parameters to hires header, but don't reslice yet.
%--------------------------------------------------------------------------

    spm_jobman('run',matlabbatch);
    
% Reset batch.
%--------------------------------------------------------------------------

    setDefaultsSPM;
    
    matlabbatch = {};
    
% Estimate registration and reslice everything to the 2mm template space.
%--------------------------------------------------------------------------

    templateDir = '/Volumes/LTD/28andMe/lpTemplate/segmentationPriors';
    
    matlabbatch{1}.spm.spatial.coreg.estwrite.ref    = cellstr(strcat([templateDir '/'], spm_select('List', templateDir, '^2mm_template.*nii$')));
    matlabbatch{1}.spm.spatial.coreg.estwrite.source = cellstr(strcat([pwd '/data.anatomical.hires/'], spm_select('List', [pwd '/data.anatomical.hires'], '^ss-ms.*nii$')));
    
    meanN3 = cellstr(strcat([pwd '/data.functional.mean/'], spm_select('List', [pwd '/data.functional.mean'], '^meanN3.*nii$')));
    rest   = cellstr(strcat([pwd '/data.functional.rest/'], spm_select('List', [pwd '/data.functional.rest'], '^ufrest4D.*nii$')));

    matlabbatch{1}.spm.spatial.coreg.estwrite.other = [meanN3; rest];
    
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep      = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol      = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm     = [7 7];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp   = 4;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap     = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask     = 0;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix   = 'r';
    
% Run the job.
%--------------------------------------------------------------------------

    spm_jobman('run', matlabbatch);
    
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