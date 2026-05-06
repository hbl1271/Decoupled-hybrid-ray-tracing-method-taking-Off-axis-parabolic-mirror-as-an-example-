%% Visualization: Scenario 2 Misalignment & Zemax Cross-Validation
%  This script simulates three critical assembly and alignment errors in the 
%  Off-Axis Parabolic (OAP) system to evaluate tolerance sensitivity:
%  1. X-Decenter (-1 mm): Lateral shift along the X-axis.
%  2. X-Tilt (-1 deg): Rotational misalignment about the X-axis.
%  3. Z-Decenter (-1 mm): Longitudinal shift along the optical axis.
%
%  Outputs include Wavefront Error (WFE) maps, Spot Diagrams, and Ray Fans 
%  for each specific error case for comparison with theoretical baselines.

clc; clear; close all;
digits(10); % Precision setting for geometric calculations

% --- Configuration & Output Setup ---
data_out_dir = fullfile(pwd, 'data');
if ~exist(data_out_dir, 'dir'), mkdir(data_out_dir); end

% Define the simulation matrix: {Identifier, Decenter_Vector [mm], Tilt_Vector [deg]}
analysis_cases = {
    'Decenter_X', [-1, 0, 0], [0, 0, 0];
    'Tilt_X',     [0, 0, 0],  [-1, 0, 0];
    'Decenter_Z', [0, 0, -1], [0, 0, 0]
};

fprintf('Starting Scenario 2: Misalignment Analysis...\n');
fprintf('==============================================\n');

for i = 1:size(analysis_cases, 1)
    case_name = analysis_cases{i,1};
    dec_val   = analysis_cases{i,2};
    tilt_val  = analysis_cases{i,3};
    
    fprintf('Running Case [%d/3]: %s\n', i, case_name);
    
    % Execute core simulation logic for the specific misalignment parameters
    results = core_misalignment_engine(dec_val, tilt_val);
    
    % Serialize results to disk for post-processing and Zemax comparison
    save_file = fullfile(data_out_dir, ['Matlab_Analysis_', case_name, '.mat']);
    save(save_file, '-struct', 'results');
    fprintf('Saved: %s\n', save_file);
end

fprintf('==============================================\n');
fprintf('Scenario 2 Complete.\n');

%% --- Core Simulation Engine ---
function out = core_misalignment_engine(decenter1, tilt1)
    % 1. OAP1 Initialization (Applying Misalignment Errors)
    D_OAP1 = 25; P_OAP1_id = [0,0,100]; F_OAP1_id = [0,0,200]; OL1 = [0,-50]; N1 = 5000;
    
    % Generate coordinates with rigid-body transformations (Decenter/Tilt)
    [Coords_OAP1_new] = generate_oap_points(D_OAP1, P_OAP1_id, F_OAP1_id, OL1, N1, decenter1, tilt1);
    % Re-calculate the actual Vertex and Focus in Global Space after transformation
    [P_OAP1_new, F_OAP1_new, ~] = fit_oap(Coords_OAP1_new, P_OAP1_id, F_OAP1_id);

    % 2. OAP2 Initialization (Ideal Stationary Mirror)
    D_OAP2 = 80; P_OAP2_id = [0,0,-100]; F_OAP2_id = [0,0,200]; OL2 = [0,-150]; N2 = 5000;
    [Coords_OAP2] = generate_oap_points(D_OAP2, P_OAP2_id, F_OAP2_id, OL2, N2);
    [R2, T2] = rigid_transform_3d(Coords_OAP2, Coords_OAP2); % Identity synchronization
    P_OAP2_new = P_OAP2_id * R2' + T2;
    F_OAP2_new = F_OAP2_id * R2' + T2;

    % 3. Source Beam Construction (High-Resolution Sampling)
    Dia_beam = 20; n_rays_axis = 256; dx = Dia_beam / n_rays_axis;
    beam_x = ((1:n_rays_axis) - n_rays_axis/2) * dx;
    beam_y = ((1:n_rays_axis) - n_rays_axis/2 - 1) * dx;
    beam_center = [0, -50];
    [Bx, By] = meshgrid(beam_x + beam_center(1), beam_y + beam_center(2));

    % Circular mask to define the pupil
    mask = ((Bx - beam_center(1)).^2 + (By - beam_center(2)).^2) <= (Dia_beam/2)^2;
    P1 = [Bx(mask), By(mask), zeros(sum(mask(:)), 1)];
    D1 = repmat([0,0,1], size(P1,1), 1);

    % Identify Chief Ray (Central Ray) for OPD referencing
    [c_grid, r_grid] = meshgrid(1:n_rays_axis, 1:n_rays_axis);
    row_v = r_grid(mask); col_v = c_grid(mask);
    mid_idx = find(row_v == n_rays_axis/2 & col_v == (n_rays_axis/2+1));

    % 4. Sequential Ray Tracing (DHRT Algorithm)
    [P2, D2] = calculate_reflection_DHRT(Coords_OAP1_new, P1, D1, P_OAP1_new, F_OAP1_new);
    [P3, D3] = calculate_reflection_DHRT(Coords_OAP2, P2, D2, P_OAP2_new, F_OAP2_new);

    % 5. Wavefront Assessment & Optical Path Difference (OPD)
    valid = ~any(isnan(P3), 2);
    P1_v = P1(valid,:); P2_v = P2(valid,:); P3_v = P3(valid,:); D3_v = D3(valid,:);
    
    L1 = vecnorm(P2_v - P1_v, 2, 2);
    L2 = vecnorm(P3_v - P2_v, 2, 2);

    % Projection to Exit Pupil (Reference distance matched to Zemax model)
    dist_exit = -2577.654; 
    D_chief = D3_v(mid_idx, :);
    P_chief_exit = P3_v(mid_idx, :) + D_chief * dist_exit;

    % Calculate segment to reference plane and final OPD
    vec_diff = P_chief_exit - P3_v;
    L3 = sum(vec_diff .* D_chief, 2) ./ sum(D3_v .* D_chief, 2);
    OPD = (L1 + L2 + L3) - (L1(mid_idx) + L2(mid_idx) + L3(mid_idx));
    Wavefront = -OPD * 1000; % Convert to micrometers

    % 6. Exit Pupil Coordinate Projection
    P_at_exit = P3_v + D3_v .* L3;
    u_ref = [1, 0, 0]; % Define local basis for the exit pupil plane
    u_vec = u_ref - sum(u_ref .* D_chief) * D_chief; u_vec = u_vec / norm(u_vec);
    v_vec = cross(D_chief, u_vec);
    P_rel = P_at_exit - P_chief_exit;
    X_exit = sum(P_rel .* u_vec, 2); Y_exit = sum(P_rel .* v_vec, 2);

    % 7. Data Normalization and Grid Interpolation
    res_grid = 256;
    x_lin = linspace(min(X_exit), max(X_exit), res_grid);
    y_lin = linspace(min(Y_exit), max(Y_exit), res_grid);
    [X_g, Y_g] = meshgrid(x_lin, y_lin);
    W_grid = griddata(X_exit, Y_exit, Wavefront, X_g, Y_g, 'natural');

    % Re-apply circular aperture mask to interpolated grid
    Dx = max(X_exit)-min(X_exit); Dy = max(Y_exit)-min(Y_exit);
    mask_p = (X_g / (Dx/2)).^2 + (Y_g / (Dy/2)).^2 <= 1.0;
    W_grid(~mask_p) = NaN;

    % 8. Consolidate Analysis for Export
    out.wavefront_data = plot_wavefront_map(x_lin, y_lin, W_grid);
    out.spot_data = analyze_and_plot_spot_diagram(Dia_beam, beam_center, ...
        Coords_OAP1_new, P_OAP1_new, F_OAP1_new, Coords_OAP2, P_OAP2_new, F_OAP2_new);
    out.rayfan_data = analyze_and_plot_ray_fan(Dia_beam, beam_center, ...
        Coords_OAP1_new, P_OAP1_new, F_OAP1_new, Coords_OAP2, P_OAP2_new, F_OAP2_new);
end