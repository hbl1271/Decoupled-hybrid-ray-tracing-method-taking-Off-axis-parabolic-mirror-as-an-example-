%% Generate a synthetic Zernike wavefront/surface error based on coefficients
function outputWavefront = generate_zernike_surface_error(D, n, a)
% GENERATE_ZERNIKE_SURFACE_ERROR Constructs a surface error map using Zernike polynomials
%
% USAGE:
%   outputWavefront = generate_zernike_surface_error(D, n, a)
%
% INPUTS:
%   D - Diameter of the aperture
%   n - Grid resolution (n x n)
%   a - Vector of Zernike coefficients [a1, a2, ..., a15]
%
% OUTPUTS:
%   outputWavefront - n x n matrix of the constructed surface
%
% DESCRIPTION:
%   This function generates a surface error map using the user-defined 
%   Zernike basis. All mathematical formulas for the 15 terms are 
%   preserved exactly as provided in the original logic.

    % Create normalized grid from -1 to 1
    coords = linspace(-1, 1, n);
    [X, Y] = meshgrid(coords, coords);
    
    % Initialize 3D array to store Zernike basis functions
    % All formulas are maintained according to the user's specific definitions
    zBasis = zeros(n, n, 15);
    
    zBasis(:,:,1)  = ones(size(X));                          % Piston
    zBasis(:,:,2)  = X;                                      % X-Tilt
    zBasis(:,:,3)  = Y;                                      % Y-Tilt
    zBasis(:,:,4)  = X.^2 + Y.^2;                            % Defocus
    zBasis(:,:,5)  = X.^2 - Y.^2;                            % Astigmatism 1
    zBasis(:,:,6)  = X.*Y;                                   % Astigmatism 2
    zBasis(:,:,7)  = (3*(X.^2 + Y.^2) - 2).*X;               % Coma X
    zBasis(:,:,8)  = (3*(X.^2 + Y.^2) - 2).*Y;               % Coma Y
    zBasis(:,:,9)  = X.*(X.^2 - 3*Y.^2);                     % Trefoil 1
    zBasis(:,:,10) = Y.*(3*X.^2 - Y.^2);                     % Trefoil 2
    zBasis(:,:,11) = (X.^2 + Y.^2).^2 - (X.^2 + Y.^2);       % Spherical
    zBasis(:,:,12) = (X.^2 - Y.^2).*(X.^2 + Y.^2);           % Higher-order 1
    zBasis(:,:,13) = X.*Y.*(X.^2 + Y.^2);                    % Higher-order 2
    zBasis(:,:,14) = X.^4 - 6*X.^2.*Y.^2 + Y.^4;             % Higher-order 3
    zBasis(:,:,15) = X.*Y.*(X.^2 - Y.^2);                    % Higher-order 4
    
    [rows, cols, numModes] = size(zBasis);
    
    % Ensure the coefficient vector has a length of 15
    coeffs = zeros(1, 15);
    lenA = min(length(a), 15);
    coeffs(1:lenA) = a(1:lenA);
    
    % Vectorized computation for efficiency: (n*n x 15) * (15 x 1)
    % Reshape the 3D basis into a 2D matrix for direct multiplication
    outputWavefront = reshape(reshape(zBasis, [], numModes) * coeffs(:), rows, cols);
end