function jacobs28andMe_calculateVDM
% Compute voxel displacement map from fieldmap scan.
%
% FORMAT jacobs28andMe_calculateVDM
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
    
% Get the phase-subtracted and short-echo magnitude images.
%--------------------------------------------------------------------------

    phaseImg = cellstr(strcat([pwd '/data.fieldmap2/'], spm_select('List', [pwd '/data.fieldmap2'], '^s.*nii$')));
    magImgs  = cellstr(strcat([pwd '/data.fieldmap1/'], spm_select('List', [pwd '/data.fieldmap1'], '^s.*nii$')));

    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.phase     = phaseImg;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.magnitude = magImgs(1);

% Specify...
%--------------------------------------------------------------------------
%   1) Short- and long-echo times for fieldmap sequence.
%   2) Blip direction (i.e. phase encoding trajectory in k-space).
%   3) Total EPI readout time (in ms).
%--------------------------------------------------------------------------

    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.et      = [4.92 7.38];
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.blipdir = 1;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.tert    = 1/16.578*1000;

% The remaining parameters can be left with SPM's defaults.
%--------------------------------------------------------------------------

    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.maskbrain       = 1;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.epifm           = 0;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.ajm             = 0;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.method   = 'Mark3D';
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.fwhm     = 10;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.pad      = 0;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.ws       = 1;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.template = cellstr(fullfile(spm('Dir'), 'toolbox/FieldMap', 'T1.nii'));
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.fwhm     = 5;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.nerode   = 2;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.ndilate  = 4;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.thresh   = 0.5;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.reg      = 0.02;
    
% Grab the first EPI so we can coregister the voxel displacement map (VDM)
% for subsequent unwarping.
%--------------------------------------------------------------------------

    epi = cellstr(strcat([pwd '/data.functional.rest/'], spm_select('List', [pwd '/data.functional.rest'], '^f.*nii$')));
    
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.session.epi   = epi(1);
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.matchvdm      = 1;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.sessname      = 'session';
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.writeunwarped = 0;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.anat          = cellstr(strcat([pwd '/data.anatomical.hires/'], spm_select('List', [pwd '/data.anatomical.hires'], '^s.*nii$')));
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.matchanat     = 0;
    
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