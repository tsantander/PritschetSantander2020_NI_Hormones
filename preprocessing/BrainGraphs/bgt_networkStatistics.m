function [networkStatistics] = bgt_networkStatistics(associationMatrix, amType, preserveEdges)
% Construct a network profile.
%
% FORMAT [networkStatistics] = bgt_networkStatistics(associationMatrix, amType, preserveEdges)
%
% REQUIRED INPUT:
%   associationMatrix
%       n x n matrix whose elements (i,j) indicate the coherence /
%       correlation between the ith and jth nodes in your network. May be
%       fully-dense or previously-thresholded based on some criterion
%       (e.g. FDR-correction).
%
%   amType
%       String specifying the type of association matrix. Enter (in
%       single quotes) either:
%           'weighted'          -   Matrix contents contain values that
%                                   indicate the strength of association
%                                   between pairs of nodes (e.g. coherence
%                                   or correlation coefficients).
%           'binary'            -   Matrix contains binary values (1,0)
%                                   that simply indicate whether an edge
%                                   exists between two nodes. Weight
%                                   information is discared.
%
%   preserveEdges
%       Indicates whether to apply an additional proportional threshold to
%       the association matrix. This value must lie between [0,1], where 1
%       indicates no additional thresholding (i.e. preserve all edges
%       present) and values < 1 indicate the percentage of edges to retain
%       (so .10 would give you the top 10% of edges, .50 would give you the
%       top 50% of edges, etc.).
%
% OUTPUT:
%   networkStatistics
%       Structure array with fields containing all diagnostic measures
%       computed below.
%__________________________________________________________________________
%
% This function will compute a series of network diagnostics given an
% association matrix (can be weighted or binary, but must be undirected).
% NOTE: depending on the size/density of your network, this can take a
% while!
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

% Initialize array to hold results, apply threshold if desired.
%--------------------------------------------------------------------------

    networkStatistics = [];
    
    if (preserveEdges < 1)
        
        associationMatrix = threshold_proportional(associationMatrix, preserveEdges);
        
    end
    
    nNodes = size(associationMatrix,2);
  
    
% Compute network statistics.
%--------------------------------------------------------------------------

    disp(['|| Estimating ' amType ' network statistics. Please wait...']);
    
    % Get nodal degree out of the way first, then compute remaining
    % statistics based on amType.
    
        networkStatistics.degrees = degrees_und(associationMatrix);

    switch amType
        
        case 'weighted'
            
            connLengths                       = weight_conversion(associationMatrix, 'lengths');
            [EBC, BC]                         = edge_betweenness_wei(connLengths);
            [D,~]                             = distance_wei(associationMatrix);
            [Eloc]                            = efficiency_wei(associationMatrix, 1);

            networkStatistics.clustering      = clustering_coef_wu(associationMatrix);
            networkStatistics.efficiencyLocal = Eloc;
            networkStatistics.transitivity    = transitivity_wu(associationMatrix);
            
            % For weighted strengths computation, determine if there are
            % any negative weights in the network (possible if weights were
            % estimated via correlation, but not an issue with coherence).
            % If so, we must calculate strengths for positive and negative
            % weights separately.
            
                if (min(associationMatrix(:)) >= 0)
                                    
                    networkStatistics.strengths = strengths_und(associationMatrix);
                    
                else
                    
                    [Spos, Sneg, vpos, vneg]                 = strengths_und_sign(associationMatrix);
                    networkStatistics.strengthsPositiveNode  = Spos;
                    networkStatistics.strengthsPositiveTotal = vpos;
                    networkStatistics.strengthsNegativeNode  = Sneg;
                    networkStatistics.strengthsNegativeTotal = vneg;
                    
                end
            
        case 'binary'
            
            [EBC, BC]                          = edge_betweenness_bin(associationMatrix);
            D                                  = distance_bin(associationMatrix);
            [Eloc]                             = efficiency_bin(associationMatrix, 1);

            networkStatistics.clustering       = clustering_coef_bu(associationMatrix);
            networkStatistics.efficiencyLocal  = Eloc;
            networkStatistics.strengths        = strengths_und(associationMatrix);
            networkStatistics.transitivity     = transitivity_bu(associationMatrix);
            
    end
    
    % Finish up with more common measures / those for which amType is
    % irrelevant.
    
        networkStatistics.betweennessNode               = BC';
        networkStatistics.betweennessEdge               = EBC;
        networkStatistics.betweennessNodeNorm           = BC / ((nNodes - 1) * (nNodes - 2));
        networkStatistics.betweennessEdgeNorm           = EBC / ((nNodes - 1) * (nNodes - 2));
        
        [lambda, Eglob, eccentricity, radius, diameter] = charpath(D);
        [density, ~, ~]                                 = density_und(associationMatrix);
        
        if (preserveEdges <= 0.50)
            [swp, ~, ~]                                 = small_world_propensity(associationMatrix);
            networkStatistics.smallWorld                = swp;
        end
        
        networkStatistics.assortativity                 = assortativity_bin(associationMatrix, 0);
        networkStatistics.charpath                      = lambda;
        networkStatistics.eccentricity                  = eccentricity;
        networkStatistics.radius                        = radius;
        networkStatistics.diameter                      = diameter;
        networkStatistics.density                       = density;
        networkStatistics.efficiencyGlobal              = Eglob;
        networkStatistics.eigenCentrality               = (eigenvector_centrality_und(associationMatrix))';
        
    % Add association matrix for reference.
    
        networkStatistics.associationMatrix             = associationMatrix;

    disp('|| Network profile constructed');
    
    
% Save output to working directory.
%--------------------------------------------------------------------------

    if (preserveEdges < 1)

        save([pwd '/networkStatistics-Top' char(num2str(preserveEdges*100)) '.mat'], 'networkStatistics');
        
    else
        
        save('networkStatistics','networkStatistics');
        
    end

end