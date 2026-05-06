%% Rigid Body Transformation (Kabsch Algorithm)
function [R, t] = rigid_transform_3d(A1, A2)
% RIGID_TRANSFORM_3D Calculates the rotation and translation between two point sets
%
% USAGE:
%   [R, t] = rigid_transform_3d(A1, A2)
%
% INPUTS:
%   A1 - Nx3 matrix of the original point cloud (one point per row)
%   A2 - Nx3 matrix of the transformed point cloud (corresponding to A1)
%
% OUTPUTS:
%   R  - 3x3 rotation matrix
%   t  - 1x3 translation vector
%
% DESCRIPTION:
%   This function finds the optimal rotation R and translation t such that 
%   (A1 * R') + t is as close as possible to A2 in a least-squares sense.
%   It uses Singular Value Decomposition (SVD) and handles reflection cases.

    % Step 1: Calculate the centroids of both point sets
    centroid1 = mean(A1, 1);  % 1x3
    centroid2 = mean(A2, 1);  % 1x3
    
    % Step 2: Center the point clouds (translate to origin)
    A1_centered = A1 - centroid1;  % Nx3
    A2_centered = A2 - centroid2;  % Nx3
    
    % Step 3: Compute the covariance matrix H
    H = A1_centered' * A2_centered;  % 3x3
    
    % Step 4: Perform Singular Value Decomposition (SVD) on H
    [U, S, V] = svd(H);
    
    % Step 5: Calculate the Rotation Matrix R
    % We must ensure R is a proper rotation matrix (det(R) = 1) 
    % and not a reflection (det(R) = -1).
    R = V * U';
    
    if det(R) < 0
        % If determinant is -1, correct the third column of V to eliminate reflection
        V(:,3) = -V(:,3);
        R = V * U';
    end
    
    % Step 6: Calculate the Translation Vector t
    % Formula derived from: centroid2 = centroid1 * R' + t
    t = centroid2 - centroid1 * R';
end