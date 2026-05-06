%% Ray-Triangle Intersection and Phong Reflected Ray Calculation
function [P2, D2, N_interp, hitMask] = calculate_reflection_Phong(Coords, P1, D1)
% CALCULATE_REFLECTION_PHONG (Smoothed Normal and Logic Corrected Version)
%
% USAGE:
%   [P2, D2, N_interp, hitMask] = calculate_reflection_Phong(Coords, P1, D1)
%
% INPUTS:
%   Coords - Nx3 matrix of surface vertex coordinates
%   P1     - Mx3 matrix of incident ray starting points
%   D1     - Mx3 matrix of incident ray directions
%
% OUTPUTS:
%   P2       - Mx3 intersection points
%   D2       - Mx3 reflected ray directions
%   N_interp - Mx3 interpolated surface normals at intersection points
%   hitMask  - Mx1 logical array indicating if a ray hit the surface
%
% LOGIC:
%   1. Computes area-weighted vertex normals from the input Coords mesh.
%   2. Performs Ray-Triangle intersection using the Möller–Trumbore algorithm.
%   3. Interpolates normals using barycentric coordinates (Phong shading logic).
%   4. Calculates reflection vectors based on the interpolated normals.

    %% 1. Data Preprocessing
    num_rays = size(P1, 1);
    num_pts = size(Coords, 1);
    
    % Normalize incident ray directions
    D1 = D1 ./ vecnorm(D1, 2, 2);
    
    % Perform Delaunay triangulation
    tri = delaunay(Coords(:,1), Coords(:,2));
    num_tris = size(tri, 1);

    %% 2. Internal Vertex Normal Calculation (Area-weighted smooth normals)
    V0_all = Coords(tri(:,1), :);
    V1_all = Coords(tri(:,2), :);
    V2_all = Coords(tri(:,3), :);
    
    % Edge vectors for area and MT algorithm
    E1 = V1_all - V0_all;
    E2 = V2_all - V0_all;
    
    % Face normals (magnitude proportional to triangle area for weighting)
    face_normals = cross(E1, E2, 2);
    
    % Accumulate face normals to vertices
    vertex_N = zeros(num_pts, 3);
    for i = 1:3
        vertex_N(:,1) = vertex_N(:,1) + accumarray(tri(:,i), face_normals(:,1), [num_pts, 1]);
        vertex_N(:,2) = vertex_N(:,2) + accumarray(tri(:,i), face_normals(:,2), [num_pts, 1]);
        vertex_N(:,3) = vertex_N(:,3) + accumarray(tri(:,i), face_normals(:,3), [num_pts, 1]);
    end
    % Normalize vertex normals
    vertex_N = vertex_N ./ vecnorm(vertex_N, 2, 2);

    %% 3. Initialize Outputs
    P2 = nan(num_rays, 3);
    D2 = nan(num_rays, 3);
    N_interp = nan(num_rays, 3);
    hitMask = false(num_rays, 1);
    epsilon = 1e-8;

    %% 4. Ray Tracing Loop (Iterate through rays against all triangles)
    for r = 1:num_rays
        orig = P1(r, :); 
        dir  = D1(r, :);
        
        % --- Step A: Filter non-parallel triangles ---
        h = cross(repmat(dir, num_tris, 1), E2, 2);
        a = dot(E1, h, 2);
        
        idx_step1 = find(abs(a) > epsilon); 
        if isempty(idx_step1), continue; end
        
        % --- Step B: Calculate and filter barycentric coordinate u ---
        f_val1 = 1 ./ a(idx_step1);
        s_val1 = orig - V0_all(idx_step1, :);
        u_val1 = f_val1 .* dot(s_val1, h(idx_step1, :), 2);
        
        check_u = (u_val1 >= 0.0) & (u_val1 <= 1.0);
        idx_step2 = idx_step1(check_u);
        if isempty(idx_step2), continue; end
        
        % --- Step C: Calculate and filter barycentric coordinate v ---
        u_val2 = u_val1(check_u);
        f_val2 = 1 ./ a(idx_step2);
        s_val2 = orig - V0_all(idx_step2, :);
        
        q = cross(s_val2, E1(idx_step2, :), 2);
        v_val2 = f_val2 .* dot(repmat(dir, length(idx_step2), 1), q, 2);
        
        check_v = (v_val2 >= 0.0) & (u_val2 + v_val2 <= 1.0);
        idx_step3 = idx_step2(check_v);
        if isempty(idx_step3), continue; end
        
        % --- Step D: Calculate and filter distance t ---
        f_val3 = 1 ./ a(idx_step3);
        q_val3 = q(check_v, :);
        t_val3 = f_val3 .* dot(E2(idx_step3, :), q_val3, 2);
        
        check_t = t_val3 > 1e-5;
        if ~any(check_t), continue; end
        
        % --- Step E: Find the closest intersection point ---
        [min_t, local_idx] = min(t_val3(check_t));
        
        % Robustly locate the winning triangle index
        valid_final_indices = idx_step3(check_t);
        best_tri_idx = valid_final_indices(local_idx);
        
        % Extract barycentric coordinates
        u_final_list = u_val2(check_v);
        u_final_valid = u_final_list(check_t);
        best_u = u_final_valid(local_idx);
        
        v_final_list = v_val2(check_v);
        v_final_valid = v_final_list(check_t);
        best_v = v_final_valid(local_idx);
        
        best_w = 1 - best_u - best_v;

        % --- Step F: Calculate physical quantities ---
        hitMask(r) = true;
        P2(r, :) = orig + min_t * dir;
        
        % Interpolate vertex normals (Phong Shading logic)
        n0 = vertex_N(tri(best_tri_idx, 1), :);
        n1 = vertex_N(tri(best_tri_idx, 2), :);
        n2 = vertex_N(tri(best_tri_idx, 3), :);
        
        n_smooth = best_w * n0 + best_u * n1 + best_v * n2;
        n_smooth = n_smooth / norm(n_smooth);
        
        % Ensure normal direction consistency
        if dot(dir, n_smooth) > 0, n_smooth = -n_smooth; end
        N_interp(r, :) = n_smooth;
        
        % Law of Reflection: R = I - 2(N.I)N
        D2(r, :) = dir - 2 * dot(dir, n_smooth) * n_smooth;
    end
    
    % Normalize all reflected ray directions
    D2 = D2 ./ vecnorm(D2, 2, 2);
end