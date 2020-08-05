function jacobs28andMe_estimateAnatomicalNoise
% Compute principal signal components in white matter + CSF voxels.
%
% FORMAT jacobs28andMe_estimateAnatomicalNoise
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
    
% Tell SPM where the hires data live, get the WM and CSF segments.
%--------------------------------------------------------------------------

    hiresPath = [pwd '/data.anatomical.hires'];
    
    whiteSeg = [hiresPath '/' spm_select('List', hiresPath, '^c2s.*nii$')];
    csfSeg   = [hiresPath '/' spm_select('List', hiresPath, '^c3s.*nii$')];
    
    cd('data.anatomical.mask');
    
% Specify inputs to ImCalc and compute mask.
%--------------------------------------------------------------------------

    matlabbatch{1}.spm.util.imcalc.input      = cellstr([whiteSeg; csfSeg]);
    matlabbatch{1}.spm.util.imcalc.output     = 'noisemask.nii';
    matlabbatch{1}.spm.util.imcalc.expression = '(i1 > 0.99) | (i2 > 0.99)';

    matlabbatch{1}.spm.util.imcalc.options.dmtx   = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask   = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 0;
    matlabbatch{1}.spm.util.imcalc.options.dtype  = 4;
    
    spm_jobman('run',matlabbatch);
    
% Register mask to mean functional scan.
%--------------------------------------------------------------------------

    meanEPI       = spm_vol(['../data.functional.mean/' spm_select('List', '../data.functional.mean', '^meanuf.*nii$')]);
    noiseMask     = spm_vol([pwd '/noisemask.nii']);
    resliceParams = struct('mean', false, 'interp', 0, 'which', 1, 'prefix', 'r');
    spm_reslice([meanEPI noiseMask], resliceParams);
    
% Navigate to where the resting-state data live.
%--------------------------------------------------------------------------

    cd('../data.functional.rest');
    
% Load in mask and 4D resting data, reshape into 2D.
%--------------------------------------------------------------------------

    rest4D = load_nii([pwd '/ufrest4D.nii.gz']);
    mask3D = load_nii('../data.anatomical.mask/rnoisemask.nii');
    
    [x,y,z,t] = size(rest4D.img);
    rest2D    = reshape(rest4D.img, x*y*z, t);
    rest2D    = double(rest2D)';
    
    mask2D    = reshape(mask3D.img, 1, numel(mask3D.img));
    
    clear rest4D mask3D
    
% Get WM/CSF voxels, normalize, and SVD.
%--------------------------------------------------------------------------

    disp('|| Extracting first 5 principal components from noise mask');

    noiseData = rest2D(:, logical(mask2D));
    zNoise    = zscore(noiseData); 
    [u,~,~]   = svd(zNoise*zNoise');
    anatNoise = u(:,1:5);
    
    disp('|| Finished component extraction. Saving...');
    
    save anatNoise anatNoise
    
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