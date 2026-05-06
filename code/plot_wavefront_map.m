%% Plot Wavefront Map and Package Visualization Data
function [plotData] = plot_wavefront_map(xLin, yLin, wGrid)
% PLOT_WAVEFRONT_MAP Visualizes the 2D wavefront error and returns metadata
%
% USAGE:
%   plotData = plot_wavefront_map(xLin, yLin, wGrid)
%
% INPUTS:
%   xLin  - 1D array of X-axis coordinates
%   yLin  - 1D array of Y-axis coordinates
%   wGrid - 2D matrix of wavefront error/sag values (typically in microns)
%
% OUTPUTS:
%   plotData - Struct containing the coordinates, grid, and calculated RMS
%
% DESCRIPTION:
%   This function renders a 2D heat map of the wavefront error. It 
%   automatically calculates the RMS value (omitting NaNs), sets up 
%   the professional color mapping (Jet), and ensures the Y-axis is 
%   oriented correctly for optical surface inspection.

    % Calculate the RMS value, omitting NaNs for masked apertures
    rmsVal = std(wGrid(:), 'omitnan');
    
    % 1. Execute Plotting
    fig = figure('Name', 'Wavefront Error'); 
    set(fig, 'Color', 'w');
    
    % Display the 2D map
    imagesc(xLin, yLin, wGrid);
    set(gca, 'Color', [0.9 0.9 0.9]); % Light gray background for NaN regions
    hold on; 
    axis image; 
    set(gca, 'YDir', 'normal'); 
    colormap('jet');
    
    % Configure colorbar and labels
    cb = colorbar; 
    ylabel(cb, '\mu m', 'FontSize', 12);
    
    % Set axis limits
    xlim([min(xLin) max(xLin)]); 
    ylim([min(yLin) max(yLin)]);
    
    % Title with dynamic RMS value display
    title(['Wavefront Error (\mum)', newline, 'RMS = ', num2str(rmsVal, '%.4f')]);
    
    % 2. Package data for storage or further processing
    plotData.x_lin = xLin;
    plotData.y_lin = yLin;
    plotData.W_grid = wGrid;
    plotData.RMS = rmsVal;
end