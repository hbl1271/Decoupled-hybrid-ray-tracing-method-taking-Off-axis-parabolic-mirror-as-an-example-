function [vFit, fFit, rmsErr] = fit_oap(surfaceCoords, vGuess, fGuess)
% FIT_OAP Fits an Off-Axis Parabolic (OAP) surface to a point cloud
%
% USAGE:
%   [vFit, fFit, rmsErr] = fit_oap(surfaceCoords, vGuess, fGuess)
%
% INPUTS:
%   surfaceCoords - Nx3 matrix of measured points
%   vGuess        - 1x3 initial guess of the vertex position
%   fGuess        - 1x3 initial guess of the focal point position
%
% OUTPUTS:
%   vFit          - 1x3 optimized vertex position
%   fFit          - 1x3 optimized focal point position
%   rmsErr        - Root-Mean-Square error of the fit (residuals)
%
% DESCRIPTION:
%   This function fixes the focal length (f) and optimizes 5 degrees of 
%   freedom: 3 translations [dx, dy, dz] and 2 rotations [rx, ry] to find 
%   the best-fit OAP geometry using non-linear least squares.

    % Calculate fixed focal length from initial guess
    fFixed = norm(fGuess - vGuess);
    
    % Optimization parameters: [dx, dy, dz, rx, ry]
    % Initialized at zero (incremental optimization)
    x0 = zeros(1, 5); 

    options = optimoptions('lsqnonlin', ...
        'Display', 'iter', ...
        'FunctionTolerance', 1e-12, ...
        'StepTolerance', 1e-12, ...
        'MaxIterations', 400);

    % Objective function: minimize the distance from points to the ideal OAP
    objFunc = @(x) calc_residuals(x, surfaceCoords, vGuess, fGuess, fFixed);

    % Perform non-linear least squares optimization
    [xOpt, resnorm] = lsqnonlin(objFunc, x0, [], [], options);

    % Recover finalized geometry
    [vFitCol, fFitCol, ~] = get_transformed_geometry(xOpt, vGuess, fGuess, fFixed);
    vFit = vFitCol';
    fFit = fFitCol';
    
    % Calculate RMS error
    rmsErr = sqrt(resnorm / size(surfaceCoords, 1));
end

%% --- Helper: Residual Calculation ---
function residuals = calc_residuals(x, coords, vInit, fInit, f)
    % 1. Update V and F based on optimization variables
    [vCurr, ~, wCurr] = get_transformed_geometry(x, vInit, fInit, f);
    
    % 2. Build local orthogonal basis [u, v, w]
    % w is the optical axis direction
    [~, minDim] = min(abs(wCurr));
    temp = [0; 0; 0]; 
    temp(minDim) = 1;
    
    u = cross(temp, wCurr); 
    u = u / norm(u);
    v = cross(wCurr, u);
    localRot = [u, v, wCurr]; 
    
    % 3. Transform point cloud into local coordinate system
    % P_local = (P_global - Vertex) * Rot_matrix
    ptsRel = coords - vCurr';
    ptsLocal = ptsRel * localRot;
    
    % 4. Calculate Parabolic Residuals
    % Local Equation: z = (x^2 + y^2) / (4f)
    zTheory = (ptsLocal(:,1).^2 + ptsLocal(:,2).^2) / (4 * f);
    residuals = ptsLocal(:,3) - zTheory;
end

%% --- Helper: Geometry Transformation ---
function [v, f, w] = get_transformed_geometry(x, vInit, fInit, focalLength)
    % Translation
    v = vInit(:) + x(1:3)';
    
    % Rotation matrices (incremental rotations)
    rx = x(4); ry = x(5);
    RotX = [1, 0, 0; 0, cos(rx), -sin(rx); 0, sin(rx), cos(rx)];
    RotY = [cos(ry), 0, sin(ry); 0, 1, 0; -sin(ry), 0, cos(ry)];
    
    % Update optical axis direction
    wInit = (fInit(:) - vInit(:)) / focalLength;
    w = (RotY * RotX * wInit); 
    
    % Update focal point position
    f = v + w * focalLength;
end