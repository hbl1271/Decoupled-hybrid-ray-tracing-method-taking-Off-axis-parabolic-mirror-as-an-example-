function [globalCoords] = generate_oap_points(D, V0, F0, OL, N_total, decenter, tilt, surface_error)
% GENERATE_OAP_POINTS Generates a point cloud for an Off-Axis Parabolic (OAP) surface
%
% USAGE:
%   Coords = generate_oap_points(D, V0, F0, OL, N_total, decenter, tilt, surface_error)
%
% INPUTS:
%   D             - Sub-aperture diameter
%   V0            - Parent parabola vertex [x, y, z]
%   F0            - Parent parabola focus [x, y, z]
%   OL            - Off-axis distance/offset [du, dv] from parent vertex to sub-aperture center
%   N_total       - Targeted total number of points
%   decenter      - Local decenter [dx, dy, dz] (default: [0,0,0])
%   tilt          - Local tilt [tx, ty, tz] in degrees (default: [0,0,0])
%   surface_error - Residual error map (2D matrix) to be interpolated
%
% OUTPUTS:
%   globalCoords  - Nx3 matrix of points in the global coordinate system

    %% 1. Initialization and Basis Construction
    if nargin < 6 || isempty(decenter), decenter = [0, 0, 0]; end
    if nargin < 7 || isempty(tilt), tilt = [0, 0, 0]; end
    if nargin < 8, surface_error = []; end

    v0_col = V0(:); f0_col = F0(:);
    optAxis = f0_col - v0_col;
    focalLength = norm(optAxis);
    w_axis = optAxis / focalLength;
    
    % Construct local-to-global rotation matrix
    temp_z = [0; 0; 1];
    u_axis = cross(temp_z, w_axis);
    if norm(u_axis) < 1e-6, u_axis = [1; 0; 0]; else u_axis = u_axis/norm(u_axis); end
    v_axis = cross(w_axis, u_axis);
    rotLocalToGlobal = [u_axis, v_axis, w_axis];
    
    %% 2. Generate Base Sub-aperture Sampling (Polar Distribution)
    R = D / 2;
    nRings = round(sqrt(N_total / pi)); 
    du_rel = 0; dv_rel = 0;
    for k = 1:nRings
        rk = R * (k / nRings); 
        nk = round(2 * pi * k); 
        theta = linspace(0, 2*pi, nk + 1);
        theta(end) = []; 
        du_rel = [du_rel, rk * cos(theta)];
        dv_rel = [dv_rel, rk * sin(theta)];
    end

    % Coordinates relative to parent parabola vertex (u, v, w)
    u_base = OL(1) + du_rel;
    v_base = OL(2) + dv_rel;
    w_base = (u_base.^2 + v_base.^2) / (4 * focalLength);
    
    %% 3. Surface Error Injection (Along Surface Normals)
    if ~isempty(surface_error)
        [rows, cols] = size(surface_error);
        % Map error grid to sub-aperture range [-R, R]
        g_u = linspace(-R, R, rows); 
        g_v = linspace(-R, R, cols);
        [G_V, G_U] = meshgrid(g_v, g_u);
        
        % Interpolate error value for each sampled point
        errorValues = interp2(G_V, G_U, surface_error, dv_rel, du_rel, 'spline');
        errorValues(isnan(errorValues)) = 0;

        % Calculate surface normals for the ideal parabola
        % n = [-dz/du, -dz/dv, 1] -> normalized
        nu = -u_base / (2 * focalLength);
        nv = -v_base / (2 * focalLength);
        nw = ones(size(u_base));
        
        norms = sqrt(nu.^2 + nv.^2 + nw.^2);
        nu = nu ./ norms; nv = nv ./ norms; nw = nw ./ norms;

        % Offset points along their local normal vectors
        u_base = u_base + errorValues .* nu;
        v_base = v_base + errorValues .* nv;
        w_base = w_base + errorValues .* nw;
    end
    
    localPoints = [u_base', v_base', w_base'];
    
    %% 4. Rigid Body Displacement (Tilt & Decenter)
    % local tilt (intrinsic rotations)
    if any(tilt ~= 0)
        rad = deg2rad(tilt);
        rx = rad(1); ry = rad(2); rz = rad(3);
        Rx = [1 0 0; 0 cos(rx) -sin(rx); 0 sin(rx) cos(rx)];
        Ry = [cos(ry) 0 sin(ry); 0 1 0; -sin(ry) 0 cos(ry)];
        Rz = [cos(rz) -sin(rz) 0; sin(rz) cos(rz) 0; 0 0 1];
        Rot = Rz * Ry * Rx;
        localPoints = localPoints * Rot';    
    end
    
    % local decenter
    localPoints = localPoints + decenter(:)';

    %% 5. Project to Global Coordinate System
    globalCoords = localPoints * rotLocalToGlobal' + v0_col';
    
end