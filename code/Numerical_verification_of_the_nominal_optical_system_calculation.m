%% Scenario 1: Algorithm Baseline & Nominal System Verification
%  This script encapsulates the simulation logic for comparing different 
%  ray-tracing methods: Möller–Trumbore (MT), Phong , and DHRT.
%  It executes iterative simulations to evaluate numerical convergence and accuracy.

clc; clear; close all;
digits(15); % Maintain high precision for coordinate calculations

% --- Path Configuration ---
% Define and create the data storage directory to ensure cross-platform compatibility
save_path = fullfile(pwd, 'data');
if ~exist(save_path, 'dir'), mkdir(save_path); end
save_full_path = fullfile(save_path, 'results.mat');

% Define the sequence of methods for the benchmark
methods = {'MT', 'Phong', 'DHRT'};
final_data = struct();

for m = 1:length(methods)
    current_method = methods{m};
    fprintf('Running simulation: %s method...\n', current_method);
    
    % --- Pre-allocate metric arrays for the current method ---
    hit_point_error = [];
    reflection_vector_error = [];
    wavefront_rms = [];

    % Iterative loop: Increases mesh density in each step (from iteration_n = 51 to 1020)
    for iteration = 1:20
        fprintf('Iteration: %d/20\n', iteration);
        iteration_n = 51 * iteration; 
        
        %% 1. OAP Surface Generation
        % Define parameters for the first OAP (Mirror 1)
        D_OAP1 = 25; P_OAP1 = [0,0,100]; F_OAP1 = [0,0,200]; OL1 = [0,0];
        [Coords_OAP1] = generate_oap_points(D_OAP1, P_OAP1, F_OAP1, OL1, iteration_n);
        
        % Define parameters for the second OAP (Mirror 2)
        D_OAP2 = 200; P_OAP2 = [0,0,-100]; F_OAP2 = [0,0,200]; OL2 = [0,0];
        [Coords_OAP2] = generate_oap_points(D_OAP2, P_OAP2, F_OAP2, OL2, iteration_n);
        
        %% 2. Laser Beam Generation (Circular Pupil Masking)
        Diameter_beam = 20;
        num_rays_axis = 51; 
        beam_range = linspace(-Diameter_beam/2, Diameter_beam/2, num_rays_axis);
        beam_mid_coords = [6,0];
        beam_x = beam_range + beam_mid_coords(1);
        beam_y = beam_range + beam_mid_coords(2);
        
        [Bx, By] = meshgrid(beam_x, beam_y);
        mask_beam = ((Bx - beam_mid_coords(1)).^2 + (By - beam_mid_coords(2)).^2) <= (Diameter_beam/2)^2;
        
        % Filter coordinates by the circular mask and define initial direction (along +Z)
        Bx = Bx(mask_beam); By = By(mask_beam); Bz = zeros(size(Bx));
        P1 = [Bx, By, Bz];
        D1 = repmat([0,0,1], size(P1,1), 1);
        
        %% 3. Ray Tracing Execution
        % Execute Sequential Ray Tracing based on the selected method's interface
        if strcmp(current_method, 'MT')
            [P2, D2] = calculate_reflection_MT(Coords_OAP1, P1, D1);
            [P3, D3] = calculate_reflection_MT(Coords_OAP2, P2, D2);
        elseif strcmp(current_method, 'Phong')
            [P2, D2] = calculate_reflection_Phong(Coords_OAP1, P1, D1);
            [P3, D3] = calculate_reflection_Phong(Coords_OAP2, P2, D2);
        elseif strcmp(current_method, 'DHRT')
            [P2, D2] = calculate_reflection_DHRT(Coords_OAP1, P1, D1, P_OAP1, F_OAP1);
            [P3, D3] = calculate_reflection_DHRT(Coords_OAP2, P2, D2, P_OAP2, F_OAP2);
        end
        
        %% 4. Wavefront Analysis Module
        % Filter out NaN values (rays that missed the surface)
        valid_mask = ~any(isnan(P1), 2) & ~any(isnan(P2), 2) & ~any(isnan(P3), 2);
        P1_valid = P1(valid_mask, :);
        P2_valid = P2(valid_mask, :);
        P3_valid = P3(valid_mask, :);
        D3_valid = D3(valid_mask, :);

        % Calculate segment lengths for Optical Path Length (OPL)
        L1_seg = vecnorm(P2_valid - P1_valid, 2, 2);
        L2_seg = vecnorm(P3_valid - P2_valid, 2, 2);

        % Project rays to a reference plane at Z = mean(P3) + 50
        D_avg = mean(D3_valid, 1);
        D_avg = D_avg / norm(D_avg);
        P_ref_center = mean(P3_valid, 1) + D_avg * 50;

        % Solve for L3 segment length via geometric intersection with the reference plane
        vec_diff = P_ref_center - P3_valid;
        numerator = sum(vec_diff .* repmat(D_avg, size(P3_valid,1), 1), 2);
        denominator = sum(D3_valid .* repmat(D_avg, size(P3_valid,1), 1), 2);
        L3_seg = numerator ./ denominator;

        % Calculate total OPL and the relative Wavefront
        OPL_total = L1_seg + L2_seg + L3_seg;
        Wavefront = mean(OPL_total) - OPL_total;

        % Remove tilt components (Piston/Tilt) using Least Squares fitting
        X_pupil = P1_valid(:, 1);
        Y_pupil = P1_valid(:, 2);
        A_fit = [X_pupil, Y_pupil, ones(size(X_pupil))];
        coeffs = A_fit \ Wavefront; 
        Wavefront_NoTilt = Wavefront - A_fit * coeffs;
        
        % Calculate RMS error (converted to micrometers)
        RMS_val = rms(Wavefront_NoTilt * 1000);

        %% 5. Ground Truth & Accuracy Evaluation
        % Obtain analytical solution for the system using matrix-based reflection
        [~, P3_ground_truth, D3_ground_truth] = get_oap_system_ground_truth_matrix...
            (P1, D1, P_OAP1, F_OAP1, P_OAP2, F_OAP2);
        
        % Compute Hit Point Error (Euclidean distance)
        hit_point_error(iteration, :) = sqrt(sum((P3 - P3_ground_truth).^2, 2));
        
        % Compute Angular Error (Arc-cosine of the dot product)
        D3_norm = D3 ./ sqrt(sum(D3.^2, 2));
        D3_gt_norm = D3_ground_truth ./ sqrt(sum(D3_ground_truth.^2, 2));
        reflection_vector_error(iteration, :) = acos(sum(D3_norm .* D3_gt_norm, 2));
        
        % Store RMS metric
        wavefront_rms(iteration) = RMS_val;
    end
    
    % Store temporary metrics into structure
    res_struct.hpe = hit_point_error;
    res_struct.rve = reflection_vector_error;
    res_struct.wr = wavefront_rms;
    
    % Organize by results1 (MT), results2 (Phong), results3 (DHRT)
    final_data.(sprintf('results%d', m)) = res_struct;
end

% --- Save Results ---
save(save_full_path, '-struct', 'final_data');
fprintf('Success: Benchmark data saved to %s\n', save_full_path);