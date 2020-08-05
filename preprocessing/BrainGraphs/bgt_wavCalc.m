function [wavScales] = bgt_wavCalc(repetitionTime, nWavelets)
% Calculation of MODWT frequency bands.
%
% FORMAT [wavScales] = bgt_wavCalc(repetitionTime, nWavelets)    
%
% REQUIRED INPUT:
%   repetitionTime
%       EPI sequence TR (entered in seconds).
%
%   nWavelets
%       Number of wavelet decompositions to perform.
%
% OUTPUT:
%   wavScales
%       nWavelet x 3 matrix containing relevant scales and their associated
%       frequency bands, ordered from highest to lowest frequency
%       harmonics.
%__________________________________________________________________________
%
% This is a helper function to determine which wavelet scales are 
% obtainable given the TR of your pulse sequence. 
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

% Initialize output.
%--------------------------------------------------------------------------

    wavScales = zeros(nWavelets,3);

    
% Compute frequency bands corresponding to the decomposition levels
% specified by nWavelets.
%--------------------------------------------------------------------------

    disp(['|| Obtaining ' num2str(nWavelets) ' wavelet scales for a signal with sampling freq ~' num2str(round(1/repetitionTime,3)) ' Hz...']);
            
    for iWav = 1:nWavelets
            
        thisScale         = [(1/(2^(iWav+1)*repetitionTime)), (1/(2^(iWav)*repetitionTime))];
        wavScales(iWav,:) = [iWav, thisScale];
        
        disp(['- Scale ' num2str(iWav) ': ~' num2str(round(thisScale(1),3)) '-' num2str(round(thisScale(2),3)) ' Hz']);
        
    end         

end