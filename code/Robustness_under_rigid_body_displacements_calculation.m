%% Scenario 2: Rigid-Body Displacements Sensitivity Analysis
%  Generates comparative plots between the DHRT method and Zemax simulations:
%  1. Wavefront Error (WFE) 3x3 Maps: Visualizes spatial phase distribution.
%  2. Zernike Coefficients Comparison: Quantifies aberrations (Tilt, Defocus, Astigmatism, etc.).
%  3. Ray Fan Distortion Plots: Evaluates transverse ray aberrations across the pupil.

clear; clc; close all;

%% 1. Configuration & Paths
data_dir = fullfile(pwd, 'data');
fig_dir  = fullfile(pwd, 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

col_titles = {'DHRT method', 'Zemax', 'Residual Error'};
row_labels = {'-1 mm decenter (X)', '-1^{\circ} tilt (X)', '-1 mm decenter (Z)'};
data_files = {
    'Matlab_Analysis_Decenter_X.mat', 'Zemax_Analysis_Decenter_X.mat';
    'Matlab_Analysis_Tilt_X.mat',    'Zemax_Analysis_Tilt_X.mat';
    'Matlab_Analysis_Decenter_Z.mat', 'Zemax_Analysis_Decenter_Z.mat'
};

% --- Professional RdBu Colormap for Residual Visualization ---
% Defines a divergent colormap to highlight positive and negative deviations
nodes = [0.00, 0.02, 0.19, 0.38; 0.20, 0.26, 0.45, 0.66; 0.40, 0.57, 0.77, 0.87; ...
         0.48, 0.88, 0.95, 0.97; 0.50, 0.97, 0.97, 0.97; 0.52, 0.96, 0.91, 0.84; ...
         0.60, 0.84, 0.70, 0.60; 0.80, 0.70, 0.24, 0.26; 1.00, 0.40, 0.00, 0.12];
x_interp = linspace(0, 1, 256);
RdBu_refined = [interp1(nodes(:,1), nodes(:,2), x_interp)', ...
                interp1(nodes(:,1), nodes(:,3), x_interp)', ...
                interp1(nodes(:,1), nodes(:,4), x_interp)'];

%% 2. Wavefront Error Comparison (3x3 Layout)
figWFE = figure('Position', [50, 50, 1300, 900], 'Color', 'w', 'Name', 'Wavefront Comparison');
t_wfe = tiledlayout(3, 3, 'TileSpacing', 'compact', 'Padding', 'loose');

panel_labels = {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)'};

for row = 1:3
    % Load Datasets for the specific misalignment case
    S_dhrt  = load(fullfile(data_dir, data_files{row, 1}));
    S_zmx   = load(fullfile(data_dir, data_files{row, 2}));
    
    W_dhrt  = flipud(S_dhrt.wavefront_data.W_grid);
    W_zmx   = S_zmx.wavefront_data.W_grid;
    W_zmx(W_zmx == 0) = NaN; % Handle Zemax invalid pixels outside the aperture

    % Resampling & Alignment (Align DHRT grid to Zemax grid resolution)
    [hz, wz] = size(W_zmx);
    if ~isequal(size(W_dhrt), [hz, wz])
        [Xz, Yz] = meshgrid(1:wz, 1:hz);
        [Xd, Yd] = meshgrid(linspace(1, wz, size(W_dhrt,2)), linspace(1, hz, size(W_dhrt,1)));
        W_dhrt_interp = interp2(Xd, Yd, W_dhrt, Xz, Yz, 'cubic');
    else
        W_dhrt_interp = W_dhrt;
    end

    % Piston removal and common mask synchronization for fair comparison
    common_mask = ~isnan(W_dhrt_interp) & ~isnan(W_zmx);
    W_dhrt_interp = W_dhrt_interp - mean(W_dhrt_interp(common_mask), 'omitnan');
    W_zmx = W_zmx - mean(W_zmx(common_mask), 'omitnan');
    
    % Calculate Residual Error and remove linear tilt components
    W_diff = W_dhrt_interp - W_zmx;
    W_diff = W_diff - mean(W_diff(common_mask), 'omitnan');
    W_diff = remove_tilt(W_diff); 

    plots = {W_dhrt_interp, W_zmx, W_diff};
    lim_val = max(abs([W_dhrt_interp(common_mask); W_zmx(common_mask)]), [], 'omitnan');

    for col = 1:3
        ax = nexttile(t_wfe);
        imagesc(plots{col}); axis image; axis off;
        
        % Manage scaling: Col 1-2 use full scale (Turbo); Col 3 uses error scale (RdBu)
        if col < 3
            colormap(ax, 'turbo');
            caxis(ax, [-lim_val, lim_val]);
        else
            colormap(ax, RdBu_refined);
            max_err = max(abs(W_diff(common_mask)), [], 'omitnan');
            caxis(ax, [-max_err, max_err]);
        end
        cb = colorbar(ax, 'eastoutside');
        title(cb, '\mu m', 'FontSize', 9);

        % Calculate and display Wavefront statistics (PV and RMS)
        v_idx = common_mask;
        rms_val = std(plots{col}(v_idx), 'omitnan');
        pv_val  = max(plots{col}(v_idx)) - min(plots{col}(v_idx));
        stats_str = sprintf('PV: %.3f \\mu m\nRMS: %.4f \\mu m', pv_val, rms_val);

        % Render Titles and Metadata
        if row == 1
            full_title = {['\fontsize{13}\bf{', col_titles{col}, '}'], '\fontsize{6} ', ['\fontsize{10}\rm{', stats_str, '}']};
            title(ax, full_title, 'Interpreter', 'tex');
        else
            title(ax, ['\fontsize{10}\rm{', stats_str, '}'], 'Interpreter', 'tex');
        end

        % Add Row descriptors (Case types)
        if col == 1
            text(-0.2, 0.5, row_labels{row}, 'Units', 'normalized', 'Rotation', 90, ...
                'HorizontalAlignment', 'center', 'FontSize', 13, 'FontWeight', 'bold');
        end
        
        % Subplot identifier labels (a-i)
        text(0, 1, panel_labels{(row-1)*3 + col}, 'Units', 'normalized', ...
            'Color','White','FontSize', 15, 'FontWeight', 'bold', 'VerticalAlignment', 'top');
    end
end
set(findall(figWFE, '-property', 'FontName'), 'FontName', 'Times New Roman');
exportgraphics(figWFE, fullfile(fig_dir, 'Fig_Wavefront_Comparison.tif'), 'Resolution', 600);

%% 3. Zernike Coefficients Decomposition (1x3 Layout)
figZernike = figure('Color', 'w', 'Position', [100, 100, 1300, 450], 'Name', 'Zernike Comparison');
t_z = tiledlayout(1, 3, 'TileSpacing', 'loose', 'Padding', 'compact');
z_panel_labels = {'(a)', '(b)', '(c)'};

for row = 1:3
    S_dhrt = load(fullfile(data_dir, data_files{row, 1}));
    S_zmx  = load(fullfile(data_dir, data_files{row, 2}));
    
    [c_dhrt, c_zmx, terms] = compute_zernike_coeffs(S_dhrt.wavefront_data.W_grid, S_zmx.wavefront_data.W_grid);

    ax = nexttile(t_z);
    Y_data = [c_dhrt(:), c_zmx(:)];
    b = bar(ax, terms, Y_data, 'grouped', 'BarWidth', 0.8);
    b(1).FaceColor = [0.2 0.45 0.7]; b(1).EdgeColor = 'none';
    b(2).FaceColor = [0.85 0.35 0.25]; b(2).EdgeColor = 'none';
    
    grid on; ax.GridLineStyle = '--'; ax.GridAlpha = 0.3;
    set(ax, 'XTick', terms, 'FontSize', 11);
    xticklabels(ax, arrayfun(@(i) sprintf('Z_{%d}', i), terms, 'UniformOutput', false));
    xlabel(ax, 'Zernike polynomial order', 'FontWeight', 'bold', 'FontSize', 12);
    if row == 1, ylabel('Zernike coefficient (\mu m)', 'FontWeight', 'bold', 'FontSize', 12); end
    title(ax, row_labels{row}, 'FontWeight', 'bold', 'FontSize', 13);
    if row == 2, legend(ax, {'DHRT method', 'Zemax'}, 'Location', 'northoutside', 'Orientation', 'horizontal'); end
    
    text(0, 1, z_panel_labels{row}, 'Units', 'normalized', ...
        'FontSize', 15, 'FontWeight', 'bold', 'VerticalAlignment', 'top');
end
set(findall(figZernike, '-property', 'FontName'), 'FontName', 'Times New Roman');
exportgraphics(figZernike, fullfile(fig_dir, 'Fig_Zernike_Comparison.tif'), 'Resolution', 600);

%% 4. Ray Fan Distortion Analysis (3x4 Layout)
figFan = figure('Color', 'w', 'Position', [100, 100, 1300, 850], 'Name', 'Ray Fan Analysis');
t_f = tiledlayout(3, 4, 'TileSpacing', 'compact', 'Padding', 'loose');
threshold = 0.01; % Threshold for dynamic scaling of small errors
fan_panel_labels = {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)', '(j)', '(k)', '(l)'};

plot_types  = {'EY_tan', 'EX_tan', 'EX_sag', 'EY_sag'};
plot_titles = {'\textbf{Tangential:} $\mathbf{E_y}$', ...
                '\textbf{Tangential:} $\mathbf{E_x}$', ...
                '\textbf{Sagittal:} $\mathbf{E_x}$', ...
                '\textbf{Sagittal:} $\mathbf{E_y}$'};
y_labels    = {'$E_y$ (mrad)', '$E_x$ (mrad)', '$E_x$ (mrad)', '$E_y$ (mrad)'};
x_labels    = {'$P_y$', '$P_y$', '$P_x$', '$P_x$'};

for row = 1:3
    S_dhrt = load(fullfile(data_dir, data_files{row, 1}));
    S_zmx  = load(fullfile(data_dir, data_files{row, 2}));
    m_rf = S_dhrt.rayfan_data;
    z_rf = S_zmx.rayfan_data;
    
    for col = 1:4
        ax = nexttile(t_f);
        var = plot_types{col};
        
        % Data mapping and Chief Ray alignment
        z_x = S_zmx.rayfan_data.(var)(:,1);
        z_y = S_zmx.rayfan_data.(var)(:,2);
        m_x = S_dhrt.rayfan_data.p_norm;
        m_y_raw = m_rf.(var); 
        [~, center_idx] = min(abs(m_x)); 
        m_y = m_y_raw - m_y_raw(center_idx); 
        
        % Visualization: Solid line for Zemax baseline, circular markers for DHRT
        p1 = plot(ax, z_x, z_y, '-', 'Color', [0.85 0.35 0.25], 'LineWidth', 1.5); hold on;
        step = 6; % Sub-sample DHRT data points for clarity
        p2 = plot(ax, m_x(1:step:end), m_y(1:step:end), 'o', 'Color', [0.15 0.45 0.75], ...
                  'MarkerSize', 4, 'LineWidth', 1.1);
        
        grid on; box on;
        set(ax, 'TickDir', 'in', 'LineWidth', 1.2, 'FontSize', 12, 'XTick', [-1, 0, 1]);
        
        % Dynamic Y-axis scaling logic to ensure readability across misalignment cases
        yl_val = max(abs(ylim(ax)));
        if yl_val < threshold
            ylim(ax, [-threshold, threshold]); 
            ax.YTick = [-threshold, 0, threshold]; 
        else
            new_limit = yl_val * 1.2;
            ylim(ax, [-new_limit, new_limit]);
            tick_val = round(new_limit, 2);
            if tick_val <= 0, tick_val = round(new_limit, 4); end
            ax.YTick = [-tick_val, 0, tick_val];
        end

        % Maintain uniform label positioning for complex subplots
        ylh = ylabel(ax, y_labels{col}, 'Interpreter', 'latex', 'FontSize', 15, 'FontWeight', 'bold');
        x_limits = xlim(ax);
        x_range = x_limits(2) - x_limits(1);
        ylh.Units = 'data'; 
        fixed_offset = -0.12 * x_range; 
        ylh.Position(1) = x_limits(1) + fixed_offset;
        ylh.HorizontalAlignment = 'center';

        if row == 1, title(ax,['\boldmath ' plot_titles{col}], 'Interpreter', 'latex', 'FontSize', 15); end
        if row == 3, xlabel(ax, ['Normalized Pupil ', x_labels{col}], 'Interpreter', 'latex'); end
        
        if col == 1
             text(-0.55, 0.5, row_labels{row}, 'Units', 'normalized', 'Rotation', 90, ...
                  'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');
        end
        if row == 1 && col == 1, lgd_items = [p1, p2]; end
        
        text(0, 0.12, fan_panel_labels{(row-1)*4 + col}, 'Units', 'normalized', ...
            'FontSize', 15, 'FontWeight', 'bold', 'VerticalAlignment', 'top');
    end
end
lgd = legend(lgd_items, {'Zemax', 'DHRT method'}, 'Orientation', 'horizontal', 'Fontsize',15);
lgd.Layout.Tile = 'north'; lgd.EdgeColor = 'none';
set(findall(figFan, '-property', 'FontName'), 'FontName', 'Times New Roman');
exportgraphics(figFan, fullfile(fig_dir, 'Fig_RayFan_Comparison.tif'), 'Resolution', 600);

%% --- Local Helper Functions ---

function [c_dhrt, c_zmx, terms] = compute_zernike_coeffs(W_dhrt_raw, W_zmx_raw)
    % Decomposes wavefront grids into Zernike coefficients using standard Fringe terms
    num_z = 9; terms = 2:num_z; % Skip Piston (Z1)
    W_dhrt = flipud(W_dhrt_raw);
    W_zmx = W_zmx_raw; W_zmx(W_zmx == 0) = NaN;
    
    [h, w] = size(W_zmx);
    [X, Y] = meshgrid(linspace(-1, 1, w), linspace(-1, 1, h));
    common_mask = ~isnan(W_zmx) & ~isnan(interp2(W_dhrt, (1:w), (1:h)', 'cubic'));
    idx = find(common_mask);
    
    % Radial and angular coordinates on the unit pupil
    r = sqrt(X(idx).^2 + Y(idx).^2); th = atan2(Y(idx), X(idx));
    
    % Basis matrix construction (Fringe Zernike 1-9)
    Z = zeros(length(idx), num_z);
    Z(:,1)=1; Z(:,2)=r.*cos(th); Z(:,3)=r.*sin(th); Z(:,4)=2*r.^2-1; 
    Z(:,5)=r.^2.*cos(2*th); Z(:,6)=r.^2.*sin(2*th); Z(:,7)=(3*r.^3-2*r).*cos(th);
    Z(:,8)=(3*r.^3-2*r).*sin(th); Z(:,9)=6*r.^4-6*r.^2+1;
    
    W_dhrt_v = interp2(W_dhrt, linspace(1,size(W_dhrt,2),w), linspace(1,size(W_dhrt,1),h)', 'cubic');
    
    % Least-squares fitting for both DHRT and Zemax data
    c_dhrt = Z \ (W_dhrt_v(idx) - mean(W_dhrt_v(idx)));
    c_zmx  = Z \ (W_zmx(idx) - mean(W_zmx(idx)));
    c_dhrt = c_dhrt(terms); c_zmx = c_zmx(terms);
end