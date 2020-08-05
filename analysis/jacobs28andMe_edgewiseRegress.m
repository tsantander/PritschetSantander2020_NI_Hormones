function [betas, tstat, pT] = jacobs28andMe_edgewiseRegress(zCoherence, zHormones)
% Edgewise regression of coherence against hormone data.
%
% FORMAT [betas, tstat, pT] = jacobs28andMe_edgewiseRegress(zCoherence, zHormones)
%
% REQUIRED INPUT:
%   zCoherence
%       nDays x nEdges matrix, with rows indexing sessions and columns
%       indexing edges in the network (normalized to unit variance)
%
%   zHormones
%       nDays x nHormones matrix, with rows indexing sessions and columns
%       indexing the hormones to be modeled (normalized to unit variance).
%
% OUTPUT:
%   betas
%       nHormones x nEdges matrix of parameter estimates (betas/slopes).
%
%   tstat
%       nHormones x nEdges matrix of t-statistics.
%
%   pT
%       nHormones x nEdges matrix of empirical p-values for t-statistics,
%       derived though 10000 iterations of nonparametric permutation
%       testing.
%__________________________________________________________________________
%
% This function regresses edgewise coherence data against hormonal
% timeseries in a time-synchronous (i.e. day-by-day) fashion. Statistical
% significance of model test-statistics is determined empirically via
% 10000 iterations of nonparametric permutation testing. 
% NOTE: rng seeds for each edge are determined by the system clock, so
% results are not perfectly replicable if this function is run multiple
% times for identical models (same coherence/hormone data).
%__________________________________________________________________________
%
% Author:
%   Tyler Santander (t.santander@psych.ucsb.edu)
%   Institute for Collaborative Biotechnologies
%   Department of Psychological & Brain Sciences
%   University of California, Santa Barbara
%   September 2019
%__________________________________________________________________________

% Preliminary setup, initialize arrays to store values.
%--------------------------------------------------------------------------

nObs   = size(zCoherence,1);
nEdge  = size(zCoherence,2);
nParam = size(zHormones,2);
betas  = zeros(nParam,nEdge);
tstat  = zeros(size(betas));
pT     = zeros(size(tstat));

% QR decomposition of design matrix (or just a vector if nHormones = 1).
% NOTE: 'perm' here is a permutation vector to ensure that abs(diag(R)) is
% decreasing during subsequent parameter assignment / standard error 
% estimation - it is NOT related to nonparametric derivation of empirical
% p-values via permutation testing.
%--------------------------------------------------------------------------

[Q,R,perm] = qr(zHormones,0);

% Compute model rank and check for deficiency (this should NOT be a problem
% if models are bivariate and/or all hormone data are standardized).
%--------------------------------------------------------------------------

if (isempty(R))
    rankX = 0;
elseif (isvector(R))
    rankX = double(abs(R(1))>0);
else
    rankX = sum(abs(diag(R)) > max(nObs,nParam)*eps(R(1)));
end
        
if (rankX < nParam)
    warning('Design matrix is rank deficient! Some regressors are collinear - interpret results with caution.');
    R    = R(1:rankX,1:rankX);
    Q    = Q(:,1:rankX);
    perm = perm(1:rankX);
end

% Loop over edges and estimate models.
%--------------------------------------------------------------------------

for iEdge = 1:nEdge
    
    % Display update every 1000 edges.
    
    if (mod(iEdge,1000) == 0), disp(['|| Edge ' num2str(iEdge) '...']); end
    
    % Get data for this edge.
    
    y = zCoherence(:,iEdge);
    
    % Estimate parameters.
    
    beta           = zeros(nParam, size(y,2));
    beta(perm,:)   = R\(Q'*y);
    betas(:,iEdge) = beta;
    
    % Obtain residuals, compute mean-squared error and standard error.

    RI         = R\eye(rankX);
    df         = max(0,nObs-rankX);
    resid      = y - zHormones*betas(:,iEdge);    
    mse        = resid'*resid/df;
    C          = (RI * RI') * mse;
    se         = zeros(nParam,1);
    se(perm,:) = sqrt(max(eps,diag(C)));
        
    % Convert betas to t-stat.
        
    tstat(:,iEdge) = betas(:,iEdge)./se;
        
    % Initialize and run NPT.
    
    rng(sum(100*clock));
    
    permB = zeros(nParam,10000);
    permT = zeros(nParam,10000);
        
    for jPerm = 1:10000
        
        % Permute data for this edge.
        
        yPerm = y(randperm(length(y)));
    
        % Estimate parameters.
    
        beta           = zeros(nParam, size(yPerm,2));
        beta(perm,:)   = R\(Q'*yPerm);
        permB(:,jPerm) = beta;
    
        % Obtain residuals, compute mean-squared error and standard error.
        
        resid      = yPerm - zHormones*permB(:,jPerm);
        mse        = resid'*resid/df;
        C          = (RI * RI') * mse;
        se         = zeros(nParam,1);
        se(perm,:) = sqrt(max(eps,diag(C)));
        
        % Convert betas to t-stat.
        
        permT(:,jPerm) = permB(:,jPerm)./se;
        
    end
    
    % Compute empirical p-values.
    
    for kParam = 1:nParam
    
        pT(kParam,iEdge) = sum(double(abs(permT(kParam,:)) >= abs(tstat(kParam,iEdge)))) / 10000;
        
    end
    
end

end