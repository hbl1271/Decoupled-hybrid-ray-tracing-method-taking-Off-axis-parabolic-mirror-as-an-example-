function rayfanData = analyze_and_plot_ray_fan(beamDiameter, beamMidCoords, coordsOAP1, pOAP1, fOAP1, coordsOAP2, pOAP2, fOAP2)
% ANALYZE_AND_PLOT_RAY_FAN Extracts tangential and sagittal ray fans and plots aberrations
%
% USAGE:
%   rayfanData = analyze_and_plot_ray_fan(beamDiameter, beamMidCoords, coordsOAP1, ...
%                pOAP1, fOAP1, coordsOAP2, pOAP2, fOAP2)
%
% INPUTS:
%   beamDiameter - Diameter of the incoming beam (mm)
%   beamMidCoords- Center coordinates of the beam [x, y]
%   coordsOAP1   - Coordinate system/Transform data for the first OAP
%   pOAP1        - Surface parameters/Vertex for the first OAP
%   fOAP1        - Focal length/Geometric parameter for the first OAP
%   coordsOAP2   - Coordinate system/Transform data for the second OAP
%   pOAP2        - Surface parameters/Vertex for the second OAP
%   fOAP2        - Focal length/Geometric parameter for the second OAP
%
% OUTPUTS:
%   rayfanData   - Structure containing normalized pupil coordinates and 
%                  angular aberrations (EX, EY) for both fans.
%
% DESCRIPTION:
%   This function performs local ray tracing for a cross-section of rays 
%   (Tangential and Sagittal). It calculates the angular deviation from 
%   the chief ray in milliradians (mrad) and generates standard ray fan plots.

    numFanRays = 51; 
    pNorm = linspace(-1, 1, numFanRays)'; 

    % --- 1. Generate Ray Source (Cross-section) ---
    % Tangential rays (Y-direction)
    p1Tan = [zeros(numFanRays, 1) + beamMidCoords(1), ...
             pNorm * (beamDiameter / 2) + beamMidCoords(2), ...
             zeros(numFanRays, 1)];
    % Sagittal rays (X-direction)
    p1Sag = [pNorm * (beamDiameter / 2) + beamMidCoords(1), ...
             zeros(numFanRays, 1) + beamMidCoords(2), ...
             zeros(numFanRays, 1)];
    
    directionIn = repmat([0, 0, 1], numFanRays, 1);

    % --- 2. Ray Tracing via DHRT (Discrete Hybrid Ray Tracing) ---
    % Trace through first OAP
    [p2Tan, d2Tan] = calculate_reflection_DHRT(coordsOAP1, p1Tan, directionIn, pOAP1, fOAP1);
    [p2Sag, d2Sag] = calculate_reflection_DHRT(coordsOAP1, p1Sag, directionIn, pOAP1, fOAP1);
    
    % Trace through second OAP
    [~, d3Tan] = calculate_reflection_DHRT(coordsOAP2, p2Tan, d2Tan, pOAP2, fOAP2);
    [~, d3Sag] = calculate_reflection_DHRT(coordsOAP2, p2Sag, d2Sag, pOAP2, fOAP2);

    % --- 3. Angular Aberration Calculation (mrad) ---
    % Normalize output directions
    d3TanNorm = d3Tan ./ vecnorm(d3Tan, 2, 2);
    d3SagNorm = d3Sag ./ vecnorm(d3Sag, 2, 2);
    
    % Reference to the chief ray (center ray of the fan)
    midIdx = ceil(numFanRays / 2);
    dChief = d3TanNorm(midIdx, :); 

    % Calculate deviations in X and Y directions
    exTan = (d3TanNorm(:, 1) - dChief(1)) * 1000;
    eyTan = (d3TanNorm(:, 2) - dChief(2)) * 1000;
    exSag = (d3SagNorm(:, 1) - dChief(1)) * 1000;
    eySag = (d3SagNorm(:, 2) - dChief(2)) * 1000;

    % --- 4. Visualization ---
    fig = figure('Name', 'Ray Fan Analysis', 'Color', 'w');
    
    % Determine plot limits for better visibility
    maxErr = max(abs([exTan; eyTan; exSag; eySag]));
    plotLim = maxErr * 1.2; 
    if plotLim == 0, plotLim = 1e-5; end 

    % Tangential Plot
    subplot(1, 2, 1);
    plot(pNorm, eyTan, 'b-', 'LineWidth', 1.5); hold on;
    plot(pNorm, exTan, 'r--', 'LineWidth', 1.5);
    grid on; set(gca, 'TickDir', 'in');
    xlim([-1, 1]); ylim([-plotLim, plotLim]);
    xlabel('Py (Normalized Pupil Y)'); ylabel('Aberration (mrad)');
    title('Tangential Fan'); legend('E_y', 'E_x', 'Location', 'best');

    % Sagittal Plot
    subplot(1, 2, 2);
    plot(pNorm, exSag, 'r-', 'LineWidth', 1.5); hold on;
    plot(pNorm, eySag, 'b--', 'LineWidth', 1.5);
    grid on; set(gca, 'TickDir', 'in');
    xlim([-1, 1]); ylim([-plotLim, plotLim]);
    xlabel('Px (Normalized Pupil X)'); ylabel('Aberration (mrad)');
    title('Sagittal Fan'); legend('E_x', 'E_y', 'Location', 'best');
    
    sgtitle('Transverse Ray Aberration');

    % --- 5. Data Packaging ---
    rayfanData.p_norm = pNorm;
    rayfanData.EX_tan = exTan;
    rayfanData.EY_tan = eyTan;
    rayfanData.EX_sag = exSag;
    rayfanData.EY_sag = eySag;
end