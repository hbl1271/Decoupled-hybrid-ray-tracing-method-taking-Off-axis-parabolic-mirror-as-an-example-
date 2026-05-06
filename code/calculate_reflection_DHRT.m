%% Decoupled Hierarchical Ray Tracing (DHRT) for OAPs
function [P2, D2] = calculate_reflection_DHRT(Coords_glo, P1_glo, D1_glo, V0_oap_glo, F0_oap_glo)
% CALCULATE_REFLECTION_DHRT Performs ray tracing by decoupling rigid-body geometry from residuals
%
% USAGE:
%   [P2, D2] = calculate_reflection_DHRT(Coords_glo, P1_glo, D1_glo, V0_oap_glo, F0_oap_glo)
%
% INPUTS:
%   Coords_glo - Nx3 point cloud of the actual OAP surface (Global)
%   P1_glo     - Mx3 incident ray positions (Global)
%   D1_glo     - Mx3 incident ray directions (Global)
%   V0_oap_glo - 3x1 or 1x3 Vertex of the OAP (Global)
%   F0_oap_glo - 3x1 or 1x3 Focus of the OAP (Global)
%
% OUTPUTS:
%   P2 - Mx3 intersection points on the real surface (Global)
%   D2 - Mx3 reflected ray directions (Global)

    %% Step 1: Local Coordinate System Setup
    F0_oap_glo = F0_oap_glo(:);
    V0_oap_glo = V0_oap_glo(:);
    V0 = V0_oap_glo';
    
    % Establish Local Basis (Local Z-axis aligned with Optical Axis)
    opt_axis = F0_oap_glo - V0_oap_glo;
    f = norm(opt_axis); 
    w_axis = opt_axis / f; 
    w_axis_temp = [0; 0; 1];
    u_axis = cross(w_axis_temp, w_axis);
    
    if norm(u_axis) < 1e-6
        u_axis = [1; 0; 0]; 
    else
        u_axis = u_axis / norm(u_axis);
    end
    v_axis = cross(w_axis, u_axis);
    
    % Transformation Matrix (Global to Local)
    Rot = [u_axis, v_axis, w_axis];
    Coords_loc = (Coords_glo - V0) * Rot;
    P1_loc = (P1_glo - V0) * Rot;
    D1_loc = D1_glo * Rot;
    V0_oap_loc = [0, 0, 0];
    F0_oap_loc = (F0_oap_glo - V0_oap_glo)' * Rot;
    
    numRays = size(P1_loc, 1);
    D1_loc = D1_loc ./ vecnorm(D1_loc, 2, 2); % Normalize incident rays
    
    % Re-calculate local parameters for clarity
    opt_axis_loc = (F0_oap_loc - V0_oap_loc)'; 
    f_loc = norm(opt_axis_loc);
    
    % Pre-allocate outputs
    P_flat = NaN(numRays, 3);
    w_res = zeros(numRays, 1);
    D2_loc = NaN(numRays, 3);
    hitFlags = false(numRays, 1);
    
    %% Step 2: Decoupling Rigid-Body Displacement and Surface Residuals
    % Ideal Parabolic Function: z = (x^2 + y^2) / 4f
    g_func = @(u,v) (u.^2 + v.^2) / (4 * f_loc);
    
    % Calculate Residuals (dw) between actual points and ideal parabola
    Coords_loc_rigid = Coords_loc;
    Coords_loc_rigid(:, 3) = g_func(Coords_loc(:, 1), Coords_loc(:, 2));
    dw = Coords_loc(:, 3) - Coords_loc_rigid(:, 3);
    
    % Triangulation of the ideal base surface
    tri = delaunay(Coords_loc_rigid(:, 1), Coords_loc_rigid(:, 2));
    numTris = size(tri, 1);
    
    % Extract triangle vertices for MT algorithm
    V1_Tri = Coords_loc_rigid(tri(:, 1), :); 
    V2_Tri = Coords_loc_rigid(tri(:, 2), :); 
    V3_Tri = Coords_loc_rigid(tri(:, 3), :);
    E1 = V2_Tri - V1_Tri;
    E2 = V3_Tri - V1_Tri;
    
    %% Step 3 & 4: Intersection and Coupled Normal Calculation
    for i = 1:numRays
        P_curr = P1_loc(i, :);
        D_curr = D1_loc(i, :);
        
        % --- MT Algorithm (Vectorized over all triangles) ---
        p_vec = cross(repmat(D_curr, numTris, 1), E2, 2); 
        det_val = sum(E1 .* p_vec, 2);
        
        keep = abs(det_val) > 1e-14;
        invDet = 1 ./ det_val;
        
        t_vec = P_curr - V1_Tri;
        k1 = sum(t_vec .* p_vec, 2) .* invDet;
        
        q_vec = cross(t_vec, E1, 2);
        k2 = sum(repmat(D_curr, numTris, 1) .* q_vec, 2) .* invDet;
        
        t = sum(E2 .* q_vec, 2) .* invDet;
        
        % Hit testing
        eps_val = 1e-9;
        hit_mask = keep & (k1 >= -eps_val) & (k2 >= -eps_val) & (k1 + k2 <= 1 + eps_val) & (t > 1e-5);

        if any(hit_mask)
            [minT, bestIdx] = min(t(hit_mask));
            hit_indices = find(hit_mask);
            TriIdx = hit_indices(bestIdx);
            hitFlags(i) = true;
            
            % 1. Intersection on the ideal surface
            P_int = P_curr + minT * D_curr;
            P_flat(i, :) = P_int;
            
            % 2. Analytic Base Slope (dz/du, dz/dv)
            k_base_u = P_int(1) / (2 * f_loc);
            k_base_v = P_int(2) / (2 * f_loc);
            
            % 3. Local Residual Slope via Least Squares fit of the triangle
            idx_v = tri(TriIdx, :);
            u_v = Coords_loc_rigid(idx_v, 1);
            v_v = Coords_loc_rigid(idx_v, 2);
            dw_v = dw(idx_v);
            
            % Solve for gradient of dw: [du, dv] * [ku; kv] = [ddw]
            A_mat = [ (u_v(2)-u_v(1)), (v_v(2)-v_v(1)); ...
                      (u_v(3)-u_v(2)), (v_v(3)-v_v(2)) ];
            B_mat = [ (dw_v(2)-dw_v(1)); ...
                      (dw_v(3)-dw_v(2)) ];
            
            k_res = A_mat \ B_mat; 
            k_res_u = k_res(1);
            k_res_v = k_res(2);
            
            % 4. Composite Normal (Analytic + Residual)
            % Normal n = unit([- (dz_base + dz_res), - (dz_base + dz_res), 1])
            n_vec = [-(k_base_u + k_res_u), -(k_base_v + k_res_v), 1];
            interpolatedNormal = n_vec / norm(n_vec); 
            
            % 5. Ensure normal orientation
            if dot(interpolatedNormal, D_curr) > 0
                interpolatedNormal = -interpolatedNormal; 
            end
            
            % 6. Reflection Calculation (Local)
            D2_loc(i, :) = D_curr - 2 * dot(D_curr, interpolatedNormal) * interpolatedNormal;
            
            % 7. Interpolate residual height at intersection (for coordinate compensation)
            bk1 = k1(TriIdx); bk2 = k2(TriIdx); bk0 = 1 - bk1 - bk2;
            w_res(i) = bk0*dw_v(1) + bk1*dw_v(2) + bk2*dw_v(3);
        end
    end
    
    % Final Coordinate Compensation (Adjust P2 from ideal surface to actual surface)
    P_rigid = P_flat;
    P_rigid(:, 3) = g_func(P_flat(:, 1), P_flat(:, 2));
    P2_loc = P_rigid;
    P2_loc(:, 3) = P2_loc(:, 3) + w_res(:);

    % Transform results back to Global Coordinate System
    P2 = P2_loc * Rot' + V0;
    D2 = D2_loc * Rot'; 
end