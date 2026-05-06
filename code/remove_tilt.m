%% Remove Linear Tilt and Piston from a Wavefront Matrix
function [wOut, tiltX, tiltY] = remove_tilt(wIn)
% REMOVE_TILT Removes linear tilt and piston terms from a wavefront matrix
%
% USAGE:
%   [wOut, tiltX, tiltY] = remove_tilt(wIn)
%
% INPUTS:
%   wIn    - Input wavefront matrix (can contain NaN values)
%
% OUTPUTS:
%   wOut   - Wavefront matrix after tilt and piston removal
%   tiltX  - Fitted slope in the X direction
%   tiltY  - Fitted slope in the Y direction
%
% DESCRIPTION:
%   This function performs a least-squares plane fit to the non-NaN 
%   elements of the input matrix. It subtracts the best-fit plane 
%   (Ax + By + C) to isolate surface residuals. Coordinates are 
%   normalized to [-1, 1] to ensure numerical stability during fitting.

    [numRows, numCols] = size(wIn);
    
    % 1. Build normalized grid coordinates (-1 to 1) 
    % Normalization prevents precision issues with large physical coordinates
    xCoords = linspace(-1, 1, numCols);
    yCoords = linspace(-1, 1, numRows);
    [X, Y] = meshgrid(xCoords, yCoords);
    
    % 2. Extract valid data points (exclude NaNs)
    validMask = ~isnan(wIn);
    wVector = wIn(validMask);
    xVector = X(validMask);
    yVector = Y(validMask);
    
    % 3. Construct the linear system: W = A*X + B*Y + C
    % Design matrix: [X, Y, 1]
    designMatrix = [xVector, yVector, ones(size(xVector))];
    
    % 4. Solve via Least Squares
    % Coefficients result in: [tiltX; tiltY; piston]
    coeffs = designMatrix \ wVector;
    
    tiltX  = coeffs(1);
    tiltY  = coeffs(2);
    piston = coeffs(3);
    
    % 5. Reconstruct the fitted tilt plane and subtract it
    wTiltPlane = coeffs(1)*X + coeffs(2)*Y + coeffs(3);
    wOut = wIn - wTiltPlane;
    
    % Maintain the original NaN mask for the output
    wOut(~validMask) = NaN;
end