function bgt_regressNuisance(timeSeries, motionParams, motionModel, addNuisance)
% Nuisance signal regression.
%
% FORMAT bgt_regressNuisance(timeSeries, nuisanceVars, nuisanceModel, addNuisance)
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
%   motionParams
%       Data can be specified here in one of two forms:
%           1) String specifying file/path information for motion
%              timeseries (i.e. translation/rotation parameters).
%              NOTE: must be a character array, not a cell array.
%           2) m x n matrix where m = number of TRs and n = number of
%              nuisance variables to regress out. 
%       NB: If you applied a temporal filter (including detrending) to the
%       BOLD data prior to this step, you MUST make sure your nuisance
%       regressors are also filtered to match the timeseries.
%
%   motionModel
%       Specify design matrix for nuisance regression. Enter (in single
%       quotes) either:
%           'basic'           -   Simple model using only the raw values
%                                 provided in motionParams.
%           'trueDerivative'  -   Uses raw values in motionParams PLUS
%                                 their temporal derivatives. Default is
%                                 to compute only the first derivative,
%                                 but this can be changed below. The
%                                 inclusion of derivatives in the design
%                                 matrix is useful because head motion
%                                 has lasting effects on the BOLD signal
%                                 beyond the TR on which it occurred
%                                 (e.g. spin history artifacts). Here,
%                                 derivatives are estimated via simple
%                                 finite differencing so that values
%                                 correspond to the same points over which
%                                 the nuisance data were originally
%                                 defined. Thus, this is a 'truer' estimate
%                                 of the gradient function.
%           'aprxDerivative'  -   Gradient estimation using backwards
%                                 differencing (i.e. for any given TR(t), 
%                                 the difference between motion at TR(t)
%                                 and TR(t-1)).
%           'AR1'             -   Time-shifted motion parameters akin to a
%                                 first order autoregressive model (i.e
%                                 predicting the effect of motion at
%                                 TR(t-1) on TR(t).
%           'fristonTrue'     -   Volterra expansion of the raw values 
%                                 provided in nuisanceVars. This includes
%                                 the parameters of the 'trueDerivative'
%                                 model PLUS their quadratic terms. These
%                                 Friston models are advantageous due to 
%                                 their ability to model known 
%                                 nonlinearities in the effect of motion on
%                                 the BOLD signal.
%           'fristonAprx'     -   Same expansion of the data as above using
%                                 the 'aprxDerivative' model.
%           'fristonAR1'      -   Expansion using the 'AR1' model. This was
%                                 the nuisance model for motion originally 
%                                 described in Friston et al. (1996).
%       NB: With a short timeseries, it's easy to over-parameterize the
%       model. We typically have tons of degrees of freedom to work with,
%       but take caution when deciding how many nuisance variables you
%       want to include.
%
% OPTIONAL INPUT:
%   addNuisance
%       Any additional nuisance variables (e.g. WM/CSF signal). Data can be
%       specified here in one of two forms:
%           1) String specifying file/path information for nuisance
%              timeseries.
%              NOTE: must be a character array, not a cell array.
%           2) m x n matrix where m = number of TRs and n = number of
%              additional nuisance variables to regress out - these are NOT
%              typically expanded like the motionParams above.
%
% OUTPUT:
%   No explicit output. If timeSeries was entered as an image volume, new
%   images will be written to the working directory; if timeSeries was
%   entered as a matrix, a new .mat file will be written to the working 
%   directory.
%__________________________________________________________________________
%
% This function performs nuisance regression of motion parameters (and
% any other signal confounds of concern).
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

% Check inputs to see if there are additional nuisance variables present.
%--------------------------------------------------------------------------

    if (nargin < 4)
        addNuisance = [];
    end


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
                t        = size(tsMatrix,1);
                
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
    
    
% Load motion parameters.
%--------------------------------------------------------------------------

    nvType = ischar(motionParams);
    
    switch nvType
        
        case 0
            
            X = motionParams;
            
        case 1
            
            X = load(motionParams);
            
    end
    
    
% Specify design matrix according to motionModel.
%--------------------------------------------------------------------------
    
    switch motionModel
        
        case 'basic'
            
            disp('|| Specifying basic model');
        
        case 'trueDerivative'
            
            disp('|| Specifying true derivative model');
            
            dx = derivative(X,1);
            X  = [X, dx];
            
        case 'aprxDerivative'
            
            disp('|| Specifying approximate derivative model');
            
            dx = [zeros(1,size(X,2)); diff(X)];
            X  = [X, dx];
            
        case 'AR1'
            
            disp('|| Specifying AR(1) model');
            
            dx = [zeros(1,size(X,2)); X(1:end-1,:)];
            X  = [X, dx];
            
        case 'fristonTrue'
            
            disp('|| Specifying Friston''s Volterra model with true derivative');
            
            dx = derivative(X,1);
            X  = [X, dx, X.^2, dx.^2];
            
        case 'fristonAprx'
            
            disp('|| Specifying Friston''s Volterra model with approximate derivative');
            
            dx = [zeros(1,size(X,2)); diff(X)];
            X  = [X, dx, X.^2, dx.^2];
            
        case 'fristonAR1'
            
            disp('|| Specifying Friston''s Volterra model with AR(1)');
            
            dx = [zeros(1,size(X,2)); X(1:end-1,:)];
            X  = [X, dx, X.^2, dx.^2];
            
        otherwise
            
            error('- ERROR: Invalid nuisance model specified -');
            
    end
    
    % Add constant.
        
        X = [ones(size(X,1),1), X];
        
% Add additional nuisance variables if present.
%--------------------------------------------------------------------------

    if (~isempty(addNuisance))

        nvType = ischar(addNuisance);
    
        switch nvType
        
            case 0
            
                X = [X, addNuisance];
            
            case 1
            
                X = [X, load(addNuisance)];
            
        end
        
    end
       
    nParam = size(X,2);
    
% Estimate nuisance models and obtain residuals.
%--------------------------------------------------------------------------

    disp(['|| Regressing out ' num2str(nParam - 1) ' nuisance covariates from ' num2str(size(tsMatrix,2)) ' features...']);
    
    startRegress = tic;

    %regressedSeries = (eye(size(tsMatrix,1)) - X*((X'*X)\X'))*tsMatrix;   %   DEPRECATED for increased numerical stability
    
    % QR decomposition of design matrix.
    
        [Q,R,perm] = qr(X,0);
    
    % Check for rank deficiency.
    
        if (isempty(R))
            p = 0;
        elseif (isvector(R))
            p = double(abs(R(1)) > 0);
        else
            p = sum(abs(diag(R)) > max(t,nParam)*eps(R(1)));
        end
        
        if (p < nParam)
            
            warning('Rank deficient nuisance matrix! :(');
            
            R    = R(1:p,1:p);
            Q    = Q(:,1:p);
            perm = perm(1:p);
        
        end
        
    % Estimate parameters, obtain residuals.
        
        betas(perm,:)   = R\(Q'*tsMatrix);
        regressedSeries = tsMatrix - X*betas;
        
    % Begin wrap-up.
    
    endRegress = toc(startRegress);
    
    disp(['|| Nuisance regression finished in ' num2str(endRegress/60) ' minutes']);
    
    timeSeries = regressedSeries; clear regressedSeries
 
    
% Write output.
%--------------------------------------------------------------------------

    switch tsType
        
        case 0
            
            save('nuisanceRegressedSeries','timeSeries');
            
        case 1
                        
            data2D              = zeros(t, x*y*z);
            data2D(:,mInd == 0) = timeSeries;
            data4D              = reshape(data2D', [x y z t]);
            
            nii.img             = data4D;
            nii.hdr             = dataHeader;
                        
            save_nii(nii, [path '/n' name ext]);
                        
    end
    
end