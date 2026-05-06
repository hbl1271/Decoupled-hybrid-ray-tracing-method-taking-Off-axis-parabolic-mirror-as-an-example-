%% Visualization: Scenario 3 Surface Deformation 
%  Generates high-fidelity comparative plots for Zernike (S1) and Modulated (S2) errors.
%  Visualizes Wavefront Error maps, 15-term Zernike histograms, and Ray Fan distortions.

clear; clc; close all;

%% 1. Configuration & Paths
data_dir = fullfile(pwd, 'data');
fig_dir  = fullfile(pwd, 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

col_titles = {'DHRT Method', 'Zemax', 'Residual Error'};
row_labels = {'{\itS}_1', '{\itS}_2'};
data_files = {
    'Matlab_Analysis_Surface_error_S1.mat', 'Zemax_Analysis_Surface_error_S1.mat';
    'Matlab_Analysis_Surface_error_S2.mat', 'Zemax_Analysis_Surface_error_S2.mat'
};

% --- Professional Divergent RdBu Colormap for Residual Visualization ---
nodes = [0.00, 0.20, 0.40, 0.48, 0.50, 0.52, 0.60, 0.80, 1.00;
         0.02, 0.26, 0.57, 0.88, 0.97, 0.96, 0.84, 0.70, 0.40;
         0.19, 0.45, 0.77, 0.95, 0.97, 0.91, 0.70, 0.24, 0.00;
         0.38, 0.66, 0.87, 0.97, 0.97, 0.84, 0.60, 0.26, 0.12]';
x_interp = linspace(0, 1, 256);
RdBu_refined = [interp1(nodes(:,1), nodes(:,2), x_interp)', ...
                interp1(nodes(:,1), nodes(:,3), x_interp)', ...
                interp1(nodes(:,1), nodes(:,4), x_interp)'];

%% 2. Wavefront Comparison (2x3 Layout)
figWFE = figure('Position', [50, 50, 1300, 750], 'Color', 'w', 'Name', 'Surface Error Wavefront Comparison');
t_wfe = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'normal');
wfe_labels = {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)'};

for row = 1:2
    % Load Comparison Datasets
    S_mat = load(fullfile(data_dir, data_files{row, 1}));
    S_zmx = load(fullfile(data_dir, data_files{row, 2}));
    
    W_mat = flipud(S_mat.wavefront_data.W_grid);
    W_zmx = S_zmx.wavefront_data.W_grid;
    W_zmx(W_zmx == 0) = NaN; % Synchronize aperture masks

    % Interpolation to match grid resolutions
    [hz, wz] = size(W_zmx);
    if ~isequal(size(W_mat), [hz, wz])
        [Xz, Yz] = meshgrid(1:wz, 1:hz);
        [Xm, Ym] = meshgrid(linspace(1, wz, size(W_mat,2)), linspace(1, hz, size(W_mat,1)));
        W_mat_interp = interp2(Xm, Ym, W_mat, Xz, Yz, 'cubic');
    else
        W_mat_interp = W_mat;
    end

    % Phase Alignment (Piston and Mean subtraction)
    common_mask = ~isnan(W_mat_interp) & ~isnan(W_zmx);
    W_mat_interp = W_mat_interp - mean(W_mat_interp(common_mask), 'omitnan');
    W_zmx = W_zmx - mean(W_zmx(common_mask), 'omitnan');
    W_diff = W_mat_interp - W_zmx;
    W_diff = W_diff - mean(W_diff(common_mask), 'omitnan');

    plots = {W_mat_interp, W_zmx, W_diff};
    lim_val = max(abs([W_mat_interp(common_mask); W_zmx(common_mask)]), [], 'omitnan');

    for col = 1:3
        ax = nexttile(t_wfe);
        imagesc(plots{col}); axis image; axis off;
        
        if col < 3
            colormap(ax, 'turbo'); % Main results in Turbo
            if ~isnan(lim_val), caxis(ax, [-lim_val, lim_val]); end
        else
            colormap(ax, RdBu_refined); % Residuals in Divergent colormap
            max_err = max(abs(W_diff(common_mask)), [], 'omitnan');
            if ~isempty(max_err), caxis(ax, [-max_err, max_err]); end
        end
        
        % Annotations and Colorbars
        cb = colorbar(ax, 'eastoutside');
        if col == 3, cb.Ruler.Exponent = 0; end
        title(cb, '\mu m', 'FontSize', 9);

        % Statistical performance metrics
        v_idx = common_mask;
        rms_v = std(plots{col}(v_idx), 'omitnan');
        pv_v = max(plots{col}(v_idx)) - min(plots{col}(v_idx));
        stats_str = sprintf('PV: %.3f \\mu m\nRMS: %.4f \\mu m', pv_v, rms_v);

        % Hierarchical Titles
        if row == 1
            full_title = {['\fontsize{13}\bf{', col_titles{col}, '}'], '\fontsize{6} ', ['\fontsize{10}\rm{', stats_str, '}']};
            title(ax, full_title, 'Interpreter', 'tex');
        else
            title(ax, ['\fontsize{10}\rm{', stats_str, '}'], 'Interpreter', 'tex');
        end

        % Row descriptors (S1 / S2)
        if col == 1
            text(-0.2, 0.5, row_labels{row}, 'Units', 'normalized', 'Rotation', 90, ...
                'HorizontalAlignment', 'center', 'FontSize', 13, 'FontWeight', 'bold');
        end
        
        % Subplot identifiers
        text(0, 1, wfe_labels{(row-1)*3 + col}, 'Units', 'normalized', ...
            'Color','White', 'FontSize', 15, 'FontWeight', 'bold', 'VerticalAlignment', 'top');
    end
end

set(findall(figWFE, '-property', 'FontName'), 'FontName', 'Times New Roman');
exportgraphics(figWFE, fullfile(fig_dir, 'Fig_Surface_Wavefront_Comparison.tif'), 'Resolution', 600);

%% 3. Zernike Decomposition (S1 - Fringe 15 Terms)
figZernike = figure('Color', 'w', 'Position', [100, 100, 900, 480], 'Name', 'Zernike Fitting Analysis');
S_mat = load(fullfile(data_dir, data_files{1, 1})); 
S_zmx = load(fullfile(data_dir, data_files{1, 2}));
[c_mat, c_zmx, terms] = compute_zernike_15(S_mat.wavefront_data.W_grid, S_zmx.wavefront_data.W_grid);

ax = axes(figZernike);
Y_data = [c_mat(:), c_zmx(:)]; 
b = bar(ax, terms, Y_data, 'grouped', 'BarWidth', 0.8);
b(1).FaceColor = [0.2 0.45 0.7]; b(1).EdgeColor = 'none';
b(2).FaceColor = [0.85 0.35 0.25]; b(2).EdgeColor = 'none';

grid on; ax.GridLineStyle = '--'; ax.GridAlpha = 0.3;
set(ax, 'XTick', terms, 'FontSize', 11);
xticklabels(ax, arrayfun(@(i) sprintf('Z_{%d}', i), terms, 'UniformOutput', false));
xlabel('Zernike polynomial order', 'FontWeight', 'bold', 'FontSize', 12); 
ylabel('Zernike coefficients (\mu m)', 'FontWeight', 'bold', 'FontSize', 12);
legend({'DHRT Method', 'Zemax'}, 'Location', 'northeast', 'Orientation', 'horizontal');

set(findall(figZernike, '-property', 'FontName'), 'FontName', 'Times New Roman');
exportgraphics(figZernike, fullfile(fig_dir, 'Fig_Surface_Zernike.tif'), 'Resolution', 600);

%% 4. Ray Fan Comparison (2x4 Layout)
figFan = figure('Color', 'w', 'Position', [100, 100, 1250, 700], 'Name', 'Surface Error Ray Fan');
t_f = tiledlayout(2, 4, 'TileSpacing', 'compact', 'Padding', 'loose');
threshold = 0.01;
fan_labels = {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)'};

for row = 1:2
    S_mat = load(fullfile(data_dir, data_files{row, 1})); 
    S_zmx = load(fullfile(data_dir, data_files{row, 2}));
    m_rf = S_mat.rayfan_data; z_rf = S_zmx.rayfan_data;
    plot_types = {'EY_tan', 'EX_tan', 'EX_sag', 'EY_sag'};
    y_labels = {'$E_y$ (mrad)', '$E_x$ (mrad)', '$E_x$ (mrad)', '$E_y$ (mrad)'};

    for col = 1:4
        ax = nexttile(t_f);
        v = plot_types{col};
        % Center Chief Ray reference
        m_x = m_rf.p_norm; m_y = m_rf.(v) - m_rf.(v)(floor(end/2)+1);
        z_x = z_rf.(v)(:,1); z_y = z_rf.(v)(:,2);

        % Comparison Plotting
        plot(ax, z_x, z_y, '-', 'Color', [0.85 0.35 0.25], 'LineWidth', 1.5); hold on;
        plot(ax, m_x(1:6:end), m_y(1:6:end), 'o', 'Color', [0.15 0.45 0.75], 'MarkerSize', 4.5);
        
        grid on; set(ax, 'TickDir', 'in', 'LineWidth', 1.2, 'FontSize', 12, 'XTick', [-1, 0, 1]);
        
        % Dynamic scaling for sensitivity
        yl = max(abs(ylim(ax)));
        if yl < threshold, yl = threshold; end
        ylim(ax, [-yl*1.2, yl*1.2]); set(ax, 'YTick', [round(-yl,2), 0, round(yl,2)]);
        
        % Axis labeling
        ylabel(ax, y_labels{col}, 'Interpreter', 'latex', 'FontSize', 14);
        if row == 2, xlabel('Normalized Pupil Coordinates', 'Interpreter', 'latex'); end
        
        if col == 1
            text(-0.55, 0.5, row_labels{row}, 'Units', 'normalized',...
                'Rotation', 90, 'HorizontalAlignment', 'center','FontSize', 14, 'FontWeight', 'bold');
        end
        text(0,1 , fan_labels{(row-1)*4 + col}, 'Units', 'normalized', ...
            'FontSize', 14, 'FontWeight', 'bold', 'VerticalAlignment', 'top');
    end
end
set(findall(figFan, '-property', 'FontName'), 'FontName', 'Times New Roman');
exportgraphics(figFan, fullfile(fig_dir, 'Fig_Surface_RayFan.tif'), 'Resolution', 600);

%% --- Helper Function ---

function [c_mat, c_zmx, terms] = compute_zernike_15(W_mat_raw, W_zmx_raw)
    % Decomposes wavefront into 15 Fringe Zernike terms
    num_z = 15; terms = 2:num_z;
    W_mat = flipud(W_mat_raw); W_zmx = W_zmx_raw; W_zmx(W_zmx == 0) = NaN;
    [h, w] = size(W_zmx); [X, Y] = meshgrid(linspace(-1, 1, w), linspace(-1, 1, h));
    
    W_mat_v = interp2(W_mat, linspace(1,size(W_mat,2),w), linspace(1,size(W_mat,1),h)', 'cubic');
    common = ~isnan(W_zmx) & ~isnan(W_mat_v);
    idx = find(common); r = sqrt(X(idx).^2 + Y(idx).^2); th = atan2(Y(idx), X(idx));
    
    % Basis Set (Fringe 1-15)
    Z = zeros(length(idx), num_z);
    Z(:,1)=1; Z(:,2)=r.*cos(th); Z(:,3)=r.*sin(th); Z(:,4)=2*r.^2-1; Z(:,5)=r.^2.*cos(2*th); 
    Z(:,6)=r.^2.*sin(2*th); Z(:,7)=(3*r.^3-2*r).*cos(th); Z(:,8)=(3*r.^3-2*r).*sin(th); 
    Z(:,9)=6*r.^4-6*r.^2+1; Z(:,10)=r.^3.*cos(3*th); Z(:,11)=r.^3.*sin(3*th);
    Z(:,12)=(4*r.^4-3*r.^2).*cos(2*th); Z(:,13)=(4*r.^4-3*r.^2).*sin(2*th);
    Z(:,14)=(10*r.^5-12*r.^3+3*r).*cos(th); Z(:,15)=(10*r.^5-12*r.^3+3*r).*sin(th);
    
    c_mat = Z \ (W_mat_v(idx) - mean(W_mat_v(idx), 'omitnan'));
    c_zmx = Z \ (W_zmx(idx) - mean(W_zmx(idx), 'omitnan'));
    c_mat = c_mat(terms); c_zmx = c_zmx(terms);
end