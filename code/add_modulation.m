function modulatedSurface = add_modulation(baseSurface, apertureD, spatialFreq, angleDeg, targetAmplitude)
% ADD_MODULATION Adds a single-frequency sinusoidal modulation to a surface profile
%
% USAGE:
%   modulatedSurface = add_modulation(baseSurface, apertureD, spatialFreq, angleDeg, targetAmplitude)
%
% INPUTS:
%   baseSurface     - Matrix (n x n) representing the initial surface profile
%   apertureD       - Clear aperture diameter (mm)
%   spatialFreq     - Spatial frequency (mm^-1), equivalent to 1/period
%   angleDeg        - Modulation orientation angle (degrees). 
%                     0 deg results in vertical fringes, 90 deg in horizontal.
%   targetAmplitude - Amplitude of the sine wave. Note: Peak-to-Valley (PV) = 2 * Amplitude.
%
% OUTPUTS:
%   modulatedSurface - The resulting surface profile with added sinusoidal error
%
% DESCRIPTION:
%   This function generates a 2D sinusoidal ripple based on a specified 
%   spatial frequency and rotation angle, then superimposes it onto the 
%   provided base surface.

    % Get grid dimensions
    [numPoints, ~] = size(baseSurface);
    
    % 1. Establish spatial coordinate system (Unit: mm)
    coords = linspace(-apertureD/2, apertureD/2, numPoints);
    [X, Y] = meshgrid(coords);
    
    % 2. Convert orientation angle to radians
    theta = deg2rad(angleDeg);
    
    % 3. Project coordinates onto the modulation direction
    % Rotate the coordinate system to align the modulation with the specified angle
    X_rot = X * cos(theta) + Y * sin(theta);
    
    % 4. Generate sinusoidal modulation
    % Formula: h = A * sin(2 * pi * f * x_rot)
    % Initial phase is set to zero for deterministic results
    modulation = targetAmplitude * sin(2 * pi * spatialFreq * X_rot);
    
    % 5. Synthesize output surface
    modulatedSurface = baseSurface + modulation;

end
