function spotData = analyze_and_plot_spot_diagram(beamDiameter, beamMidCoords, coordsOAP1, pOAP1, fOAP1, coordsOAP2, pOAP2, fOAP2)
% ANALYZE_AND_PLOT_SPOT_DIAGRAM Generates a geometric spot diagram using hexapolar sampling
%
% USAGE:
%   spotData = analyze_and_plot_spot_diagram(beamDiameter, beamMidCoords, ...
%              coordsOAP1, pOAP1, fOAP1, coordsOAP2, pOAP2, fOAP2)
%
% INPUTS:
%   beamDiameter - Diameter of the incoming beam (mm)
%   beamMidCoords- Center coordinates of the beam [x, y]
%   coordsOAP1   - Coordinate system/Transform data for the first OAP
%   pOAP1        - Surface parameters for the first OAP
%   fOAP1        - Geometric parameters for the first OAP
%   coordsOAP2   - Coordinate system/Transform data for the second OAP
%   pOAP2        - Surface parameters for the second OAP
%   fOAP2        - Geometric parameters for the second OAP
%
% OUTPUTS:
%   spotData     - Structure containing validated spot coordinates (SX, SY)
%                  and calculated RMS/GEO radii.
%
% DESCRIPTION:
%   The function uses hexapolar (ring-based) ray sampling to simulate a 
%   uniform beam cross-section. It traces rays through the system and 
%   plots the angular distribution (spot diagram) in milliradians.

    % --- 1. Generate Hexapolar Ray Source ---
    numRings = 10; 
    hexX = 0; hexY = 0;
    for r = 1:numRings
        numPts = r * 6;
        rho = (r / numRings) * (beamDiameter / 2); 
        theta = linspace(0, 2*pi, numPts + 1);
        theta(end) = []; % Remove redundant end point
        hexX = [hexX, cos(theta) * rho];
        hexY = [hexY, sin(theta) * rho];
    end
    
    p1Hex = [hexX' + beamMidCoords(1), hexY' + beamMidCoords(2), zeros(length(hexX), 1)];
    directionIn = repmat([0, 0, 1], size(p1Hex, 1), 1);

    % --- 2. Ray Tracing via DHRT ---
    [p2Hex, d2Hex] = calculate_reflection_DHRT(coordsOAP1, p1Hex, directionIn, pOAP1, fOAP1);
    [~, d3Hex] = calculate_reflection_DHRT(coordsOAP2, p2Hex, d2Hex, pOAP2, fOAP2);

    % --- 3. Angular Aberration Calculation (mrad) ---
    d3HexNorm = d3Hex ./ vecnorm(d3Hex, 2, 2);
    dChief = d3HexNorm(1, :); % First point assumed as chief ray
    
    sxMrad = (d3HexNorm(:, 1) - dChief(1)) * 1000;
    syMrad = (d3HexNorm(:, 2) - dChief(2)) * 1000;

    % Filter valid intersections
    validIdx = ~isnan(sxMrad) & ~isnan(syMrad);
    sxFinal = sxMrad(validIdx);
    syFinal = syMrad(validIdx);

    % --- 4. Statistical Analysis ---
    centroidX = mean(sxFinal);
    centroidY = mean(syFinal);
    radialDist = sqrt((sxFinal - centroidX).^2 + (syFinal - centroidY).^2);
    rmsRadius = sqrt(mean(radialDist.^2));
    geoRadius = max(radialDist);

    % --- 5. Visualization ---
    fig = figure('Name', 'Spot Diagram', 'Color', 'w');
    plot(sxFinal, syFinal, 'b+', 'MarkerSize', 5, 'LineWidth', 0.6); 
    axis image; grid on;
    
    % Set scale limits (0.2 mrad total window)
    limVal = 0.1;
    xlim([-limVal, limVal]); 
    ylim([-limVal, limVal]);
    
    % Style the axes as a scale indicator
    % NOTE: XTick and YTick should be adjusted dynamically based on spot size.
    set(gca, 'GridColor', [0.8 0.8 0.8], ...
        'XTick', [-0.5, 0.5], ... % Hide ticks outside range
        'YTick', [-0.5, 0.5], ...
        'XTickLabel', {}, ...      
        'YTickLabel', {}, ...      
        'TickDir', 'in',...
        'LineWidth', 1.5,...
        'FontSize', 12,...
        'FontName', 'Times New Roman');
    
    % Display statistics box
    infoStr = sprintf('RMS radius: %.4f mrad\nGEO radius: %.4f mrad', rmsRadius, geoRadius);
    text(0.095, 0.095, infoStr, ...
        'FontSize', 11, 'FontName', 'Times New Roman', 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
        'BackgroundColor', [1 1 1 0.8], 'EdgeColor', [0.5 0.5 0.5]);
    
    % Minimalist axis styling (Clean scale-bar look)
    box off;
    set(gca, 'XColor', 'none');       
    set(gca, 'YColor', 'k');          
    ylabel('0.2 mrad Scale', 'FontWeight', 'Bold'); 
    
    % --- 6. Data Packaging ---
    spotData.SX_final = sxFinal;
    spotData.SY_final = syFinal;
    spotData.RMS_radius = rmsRadius;
    spotData.GEO_radius = geoRadius;
end