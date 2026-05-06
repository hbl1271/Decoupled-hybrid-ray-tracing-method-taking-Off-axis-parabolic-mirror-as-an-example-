%% Scenario 3: Complex Surface Deformations Analysis
%  This script evaluates the system's sensitivity to high-fidelity surface 
%  topography errors on the primary OAP mirror. It performs:
%  1. Generation of Zernike-based (S1) and Modulated (S2) surface errors.
%  2. Export of surface data to Zemax-compatible Grid Sag formats.
%  3. Sequential ray tracing using the DHRT algorithm with perturbed normals.
%  4. Comparative visualization of PV/RMS statistics.

clc; clear; close all;
digits(10);

% --- 0. Path & Configuration ---
data_dir = fullfile(pwd, 'data');
zemax_dir = fullfile(pwd, 'zemax_data');
fig_dir  = fullfile(pwd, 'figures');

if ~exist(data_dir, 'dir'), mkdir(data_dir); end
if ~exist(zemax_dir, 'dir'), mkdir(zemax_dir); end
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

% Mirror and Sampling Parameters
D_aperture = 25;        % Pupil diameter (mm)
N_grid = 1000;          % Surface sampling density
decenter_origin = [0, -50]; % Off-axis decenter coordinates
lin_coords = linspace(-D_aperture/2, D_aperture/2, N_grid);
[X, Y] = meshgrid(lin_coords, lin_coords);
mask = (sqrt(X.^2 + Y.^2) <= (D_aperture / 2));

% Tight bounding box for visualization
[rows, cols] = find(mask);
row_range = min(rows):max(rows);
col_range = min(cols):max(cols);

%% --- 1. Surface Generation (Low-Frequency vs. Mid-Frequency) ---
fprintf('Step 1: Generating complex surface error profiles...\n');

% Case S1: Low-order aberrations modeled by Zernike coefficients
% Represents typical manufacturing figure errors (Trefoil, Astigmatism, Coma)
z_coefs = [0, 0, 0, 0.4, 0.4, 0.4, 0.4, 0.4, 0.5, 0.5, 0.6, 0.4, 0, 0.4, 0.4];
S1_raw = generate_zernike_surface_error(D_aperture, N_grid, z_coefs) * 1e-4;

% Case S2: Mid-frequency modulation (Sinusoidal ripple)
% Represents tool marks or periodic errors from diamond turning
S2_error_base = zeros(N_grid, N_grid);
S2_raw = add_modulation(S2_error_base, D_aperture, 0.2, 45, 1e-4);

% Data Preparation for Plotting
crop_data1 = S1_raw(row_range, col_range) * 1000; % Convert to microns
crop_data2 = S2_raw(row_range, col_range) * 1000;
crop_mask  = mask(row_range, col_range);

crop_data1(~crop_mask) = NaN;
crop_data2(~crop_mask) = NaN;
data_list = {crop_data1, crop_data2};

% --- 2. Visualization of Surface Topography ---
fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [5, 5, 22, 11]);
t = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'loose');
labels = {'(a)', '(b)'};

for col = 1:2
    ax = nexttile(t);
    imagesc(data_list{col});
    axis image; axis off;
    colormap(ax, 'turbo');
    
    current_data = data_list{col}(crop_mask);
    if ~isempty(current_data)
        local_lim = max(abs(current_data), [], 'omitnan');
        if local_lim > 0, caxis(ax, [-local_lim, local_lim]); end
    end

    % Statistical Summary (PV and RMS)
    pv = max(current_data) - min(current_data);
    rms = std(current_data);
    stats_str = sprintf('PV: %.3f \\mu m\nRMS: %.4f \\mu m', pv, rms);

    % Professional Figure Labeling
    title(ax, ['\fontsize{10}\rm{', stats_str, '}'], 'Interpreter', 'tex');
    cb = colorbar(ax, 'eastoutside');
    cb.Ruler.Exponent = 0; 
    title(cb, '\mu m', 'FontSize', 9);
    
    text(ax, 0, 1, labels{col}, 'Units', 'normalized', 'Color', 'w', ...
        'FontSize', 15, 'FontWeight', 'bold', 'VerticalAlignment', 'top');
end

set(findall(fig, '-property', 'FontName'), 'FontName', 'Times New Roman');
exportgraphics(fig, fullfile(fig_dir, 'Fig_Surface_Analysis.tif'), 'Resolution', 600);

%% --- 3. Data Export & Ray Tracing Simulation ---

% Export to Zemax Grid Sag format for cross-validation
matrix2zemax(S1_raw, D_aperture, decenter_origin(1), decenter_origin(2), fullfile(zemax_dir, 'S1_surface.dat'));
matrix2zemax(S2_raw, D_aperture, decenter_origin(1), decenter_origin(2), fullfile(zemax_dir, 'S2_surface.dat'));

% Restore matrix orientation for the DHRT engine (aligning with coordinate convention)
S1_mat = fliplr(S1_raw); 
S2_mat = fliplr(S2_raw);

surfaces = {S1_mat, S2_mat};
case_names = {'Surface_error_S1', 'Surface_error_S2'};

for i = 1:2
    fprintf('   Processing Simulation Case %d: %s\n', i, case_names{i});
    results = core_surface_analysis_engine(surfaces{i});
    save_path = fullfile(data_dir, ['Matlab_Analysis_', case_names{i}, '.mat']);
    save(save_path, '-struct', 'results');
end

%% --- Core Simulation Engine (Surface Error Integration) ---
function out = core_surface_analysis_engine(input_surface_error)
    % 1. OAP1 Setup: Integrated Surface Sag Error
    D1=25; P1_id=[0,0,100]; F1_id=[0,0,200]; OL1=[0,-50]; N1=5000;
    [Coords_OAP1_new] = generate_oap_points(D1, P1_id, F1_id, OL1, N1, [0,0,0], [0,0,0], input_surface_error);
    [P1_fit, F1_fit, ~] = fit_oap(Coords_OAP1_new, P1_id, F1_id);

    % 2. OAP2 Setup: Nominal System Geometry
    D2=200; P2_id=[0,0,-100]; F2_id=[0,0,200]; OL2=[0,-150]; N2=5000;
    [Coords_OAP2] = generate_oap_points(D2, P2_id, F2_id, OL2, N2);
    [R2, T2] = rigid_transform_3d(Coords_OAP2, Coords_OAP2);
    P2_fit = P2_id * R2' + T2; F2_fit = F2_id * R2' + T2;

    % 3. Source Beam Initialization
    dia_beam = 20; n_axis = 256; dx = dia_beam / n_axis;
    b_x = ((1:n_axis) - n_axis/2) * dx; b_y = ((1:n_axis) - n_axis/2 - 1) * dx;
    b_center = [0, -50];
    [Bx, By] = meshgrid(b_x + b_center(1), b_y + b_center(2));
    mask = ((Bx - b_center(1)).^2 + (By - b_center(2)).^2) <= (dia_beam/2)^2;
    P_src = [Bx(mask), By(mask), zeros(sum(mask(:)), 1)];
    D_src = repmat([0,0,1], size(P_src,1), 1);

    % Reference Chief Ray Index
    [c_g, r_g] = meshgrid(1:n_axis, 1:n_axis);
    row_v = r_g(mask); col_v = c_g(mask);
    mid_ptr = find(row_v == n_axis/2 & col_v == (n_axis/2+1));

    % 4. Sequential Simulation (DHRT Method)
    [P2_t, D2_t] = calculate_reflection_DHRT(Coords_OAP1_new, P_src, D_src, P1_fit, F1_fit);
    [P3_t, D3_t] = calculate_reflection_DHRT(Coords_OAP2, P2_t, D2_t, P2_fit, F2_fit);

    % 5. Wavefront & OPD Calculation
    valid = ~any(isnan(P3_t), 2);
    P3_v = P3_t(valid,:); D3_v = D3_t(valid,:); P2_v = P2_t(valid,:); P_src_v = P_src(valid,:);
    L1 = vecnorm(P2_v - P_src_v, 2, 2);
    L2 = vecnorm(P3_v - P2_v, 2, 2);

    dist_exit = -2074.926; % Reference plane for exit pupil
    D_chief = D3_v(mid_ptr, :);
    P_chief_exit = P3_v(mid_ptr, :) + D_chief * dist_exit;

    vec_diff = P_chief_exit - P3_v;
    L3 = sum(vec_diff .* D_chief, 2) ./ sum(D3_v .* D_chief, 2);
    OPD = (L1 + L2 + L3) - (L1(mid_ptr) + L2(mid_ptr) + L3(mid_ptr));
    Wavefront = -OPD * 1000; 

    % 6. Coordinate Mapping to Pupil Plane
    P_at_exit = P3_v + D3_v .* L3;
    u_ref = [1, 0, 0];
    u_vec = u_ref - sum(u_ref .* D_chief) * D_chief; u_vec = u_vec / norm(u_vec);
    v_vec = cross(D_chief, u_vec);
    P_rel = P_at_exit - P_chief_exit;
    X_exit = sum(P_rel .* u_vec, 2); Y_exit = sum(P_rel .* v_vec, 2);

    % 7. Final Resampling and Tilt Removal
    res_grid = 256;
    x_lin = linspace(min(X_exit), max(X_exit), res_grid);
    y_lin = linspace(min(Y_exit), max(Y_exit), res_grid);
    [X_g, Y_g] = meshgrid(x_lin, y_lin);
    W_grid = griddata(X_exit, Y_exit, Wavefront, X_g, Y_g, 'natural');

    Dx = max(X_exit)-min(X_exit); Dy = max(Y_exit)-min(Y_exit);
    mask_p = (X_g / (Dx/2)).^2 + (Y_g / (Dy/2)).^2 <= 1.0;
    W_grid(~mask_p) = NaN;
    W_grid_no_tilt = remove_tilt(W_grid);
    
    % 8. Package Results for Visualization
    out.wavefront_data.W_grid = W_grid_no_tilt;
    out.wavefront_data.x_lin = x_lin;
    out.wavefront_data.y_lin = y_lin;
    out.spot_data = analyze_and_plot_spot_diagram(dia_beam, b_center, ...
        Coords_OAP1_new, P1_fit, F1_fit, Coords_OAP2, P2_fit, F2_fit);
    out.rayfan_data = analyze_and_plot_ray_fan(dia_beam, b_center, ...
        Coords_OAP1_new, P1_fit, F1_fit, Coords_OAP2, P2_fit, F2_fit);
end