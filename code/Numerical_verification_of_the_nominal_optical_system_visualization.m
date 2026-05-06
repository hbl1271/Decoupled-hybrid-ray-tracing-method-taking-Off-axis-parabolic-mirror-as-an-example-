%% Visualization: Scenario 1 Baseline Results
%  Consolidates results from benchmarking simulations and generates comparative plots.
%  This script is optimized for publication-quality output with semi-logarithmic 
%  scaling to visualize convergence characteristics across several orders of magnitude.

clear; clc; close all;

%% 1. Data Loading and Unit Initialization
data_folder = fullfile(pwd, 'data');
file_name   = 'results.mat';
load_path   = fullfile(data_folder, file_name);

if ~exist(load_path, 'file')
    error('Data file not found. Please run the simulation script first.');
end

data = load(load_path);

% Unit Conversion: Preserving original logic for directional error
% Converts specific internal units to milliradians (mrad)
data.results1.rve = data.results1.rve / 360 * 2 * pi * 1000;
data.results2.rve = data.results2.rve / 360 * 2 * pi * 1000;
data.results3.rve = data.results3.rve / 360 * 2 * pi * 1000;

% X-axis represents the number of sampling points (mesh density)
x_axis = (1:20) * 51; 
colors = [242, 53, 87; 253, 185, 107; 34, 178, 218] / 255; % Publication-ready color palette
m_idx  = 2:2:length(x_axis); % Marker indices for visual clarity (every 2nd point)

%% 2. Figure Configuration
figCombined = figure(10);
% Set physical dimensions (centimeters) for consistent export scaling
set(figCombined, 'Units', 'centimeters', 'Position', [2, 2, 28, 10]); 

tlo = tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'loose');

%% --- Subplot (a): Mean Intersection Error ---
% Evaluates how accurately each method finds the ray-surface intersection point
nexttile;
p1 = semilogy(x_axis, mean(data.results1.hpe, 2, 'omitnan'), '-o', 'Color', colors(1,:), 'LineWidth', 1.5, 'MarkerSize', 5, 'MarkerIndices', m_idx, 'MarkerFaceColor', 'w'); hold on;
p2 = semilogy(x_axis, mean(data.results2.hpe, 2, 'omitnan'), '--^', 'Color', colors(2,:), 'LineWidth', 1.5, 'MarkerSize', 5, 'MarkerIndices', m_idx, 'MarkerFaceColor', 'w');
p3 = semilogy(x_axis, mean(data.results3.hpe, 2, 'omitnan'), ':d', 'Color', colors(3,:), 'LineWidth', 1.5, 'MarkerSize', 5, 'MarkerIndices', m_idx, 'MarkerFaceColor', 'w');

ylabel('Mean intersection error (mm)', 'FontWeight', 'bold');
title('(a)', 'Units', 'normalized', 'FontSize', 14, 'FontName', 'Times New Roman', 'HorizontalAlignment', 'left');
set_std_axis(gca);
set_dense_ticks(gca);

%% --- Subplot (b): Directional (Angular) Error ---
% Evaluates the accuracy of the calculated reflection vectors compared to ground truth
nexttile;
semilogy(x_axis, mean(data.results1.rve, 2, 'omitnan'), '-o', 'Color', colors(1,:), 'LineWidth', 1.5, 'MarkerSize', 5, 'MarkerIndices', m_idx, 'MarkerFaceColor', 'w'); hold on;
semilogy(x_axis, mean(data.results2.rve, 2, 'omitnan'), '--^', 'Color', colors(2,:), 'LineWidth', 1.5, 'MarkerSize', 5, 'MarkerIndices', m_idx, 'MarkerFaceColor', 'w');
semilogy(x_axis, mean(data.results3.rve, 2, 'omitnan'), ':d', 'Color', colors(3,:), 'LineWidth', 1.5, 'MarkerSize', 5, 'MarkerIndices', m_idx, 'MarkerFaceColor', 'w');

ylabel('Mean directional error (mrad)', 'FontWeight', 'bold');
title('(b)', 'Units', 'normalized', 'FontSize', 14, 'FontName', 'Times New Roman', 'HorizontalAlignment', 'left');
set_std_axis(gca);
set_dense_ticks(gca); 

%% --- Subplot (c): Wavefront RMS Error ---
% Evaluates the system-wide optical performance and phase consistency
nexttile;
semilogy(x_axis, data.results1.wr(:), '-o', 'Color', colors(1,:), 'LineWidth', 1.5, 'MarkerSize', 5, 'MarkerIndices', m_idx, 'MarkerFaceColor', 'w'); hold on;
semilogy(x_axis, data.results2.wr(:), '--^', 'Color', colors(2,:), 'LineWidth', 1.5, 'MarkerSize', 5, 'MarkerIndices', m_idx, 'MarkerFaceColor', 'w');
semilogy(x_axis, data.results3.wr(:), ':d', 'Color', colors(3,:), 'LineWidth', 1.5, 'MarkerSize', 5, 'MarkerIndices', m_idx, 'MarkerFaceColor', 'w');

ylabel('RMS wavefront error (\mum)', 'FontWeight', 'bold');
title('(c)', 'Units', 'normalized', 'FontSize', 14, 'FontName', 'Times New Roman', 'HorizontalAlignment', 'left');
set_std_axis(gca);
set_dense_ticks(gca);

%% 3. Global Labels & Legend
lgd = legend([p1, p2, p3], {'Möller–Trumbore method', 'Phong method', 'DHRT method'});
lgd.Layout.Tile = 'north'; 
lgd.Orientation = 'horizontal';
set(lgd, 'Box', 'off', 'FontSize', 12, 'FontName', 'Times New Roman');

xlabel(tlo, 'Number of sampling points on OAP', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Times New Roman');

%% 4. Figure Export
% Saves high-resolution TIFF for inclusion in academic manuscripts
fig_output_path = fullfile(pwd, 'figures');
if ~exist(fig_output_path, 'dir'), mkdir(fig_output_path); end
exportgraphics(figCombined, fullfile(fig_output_path, 'Fig2_Optimized.tif'), 'Resolution', 600);

%% --- Helper Functions for Aesthetic Consistency ---

function set_std_axis(ax)
    % Configure standard plot parameters for scientific aesthetics
    grid(ax, 'off');
    ax.YGrid = 'on';
    ax.GridLineStyle = ':';
    ax.GridAlpha = 0.3;
    ax.TickDir = 'out';
    ax.LineWidth = 0.8;
    ax.FontName = 'Times New Roman';
    ax.FontSize = 12;
    box(ax, 'off');
end

function set_dense_ticks(ax)
    % Custom logarithmic tick management for wide dynamic range data
    y_limits = ylim(ax);
    y_min = y_limits(1);
    y_max = y_limits(2);    
    upper_exp = ceil(log10(y_max));
    lower_exp = floor(log10(y_min));
    % Generate ticks at every second power of 10 for cleaner display
    tick_exponents = upper_exp:-2:(lower_exp-1);
    tick_values = 10.^tick_exponents;
    ax.YTick = sort(tick_values);
    ylim(ax, [y_min, 10^upper_exp]);
    set(ax, 'YScale', 'log'); 
end