function [intersectPoints, reflectDirs] = calculate_reflection_MT(surfaceCoords, rayOrigins, rayDirections)
% CALCULATE_REFLECTION_MT Ray tracing reflection using the Möller-Trumbore algorithm
%
% USAGE:
%   [P2, D2] = calculate_reflection_MT(surfaceCoords, rayOrigins, rayDirections)
%
% INPUTS:
%   surfaceCoords - Nx3 matrix of surface point coordinates
%   rayOrigins    - Mx3 matrix of incident ray starting positions
%   rayDirections - Mx3 matrix of incident ray directions
%
% OUTPUTS:
%   intersectPoints - Mx3 matrix of intersection points on the mesh
%   reflectDirs     - Mx3 matrix of reflected ray directions
%
% DESCRIPTION:
%   This function implements the standard Möller-Trumbore ray-triangle 
%   intersection algorithm. Unlike the hybrid DHRT method, this version 
%   calculates reflection based on the discrete geometric face normal of 
%   the hit triangle (Edge1 x Edge2) without surface gradient interpolation.

    %% 1. Pre-processing
    numRays = size(rayOrigins, 1);
    epsilon = 1e-8; % Tolerance for parallel ray-plane test
    
    % Normalize incident ray directions
    rayDirections = rayDirections ./ vecnorm(rayDirections, 2, 2);
    
    % Perform Delaunay triangulation for mesh generation
    triIndices = delaunay(surfaceCoords(:, 1), surfaceCoords(:, 2));    
    numTris = size(triIndices, 1);
    
    % Extract triangle vertices
    v0 = surfaceCoords(triIndices(:, 1), :);
    v1 = surfaceCoords(triIndices(:, 2), :);
    v2 = surfaceCoords(triIndices(:, 3), :);
    
    % Compute edge vectors
    e1 = v1 - v0;
    e2 = v2 - v0;
    
    % Calculate discrete face normals (Normalized Edge1 x Edge2)
    faceNormals = cross(e1, e2, 2);
    faceNormals = faceNormals ./ vecnorm(faceNormals, 2, 2); 
    
    %% 2. Initialization
    intersectPoints = nan(numRays, 3);
    reflectDirs = nan(numRays, 3);
    
    %% 3. Ray Tracing Loop
    for r = 1:numRays
        orig = rayOrigins(r, :);
        dir  = rayDirections(r, :);
        
        % --- Moller-Trumbore Ray-Triangle Intersection Test ---
        h = cross(repmat(dir, numTris, 1), e2, 2);
        a = dot(e1, h, 2);
        
        % Filter out rays parallel to triangles
        validTrisIdx = find(abs(a) > epsilon);
        if isempty(validTrisIdx), continue; end
        
        f = 1 ./ a(validTrisIdx);
        s = orig - v0(validTrisIdx, :);
        u = f .* dot(s, h(validTrisIdx, :), 2);
        
        % Barycentric coordinate u check
        checkU = (u >= 0.0) & (u <= 1.0);
        if ~any(checkU), continue; end
        
        % Narrow down to triangles passing u-test
        idxStepU = validTrisIdx(checkU);
        f = f(checkU);
        s = s(checkU, :);
        u = u(checkU);
        
        q = cross(s, e1(idxStepU, :), 2);
        v = f .* dot(repmat(dir, length(f), 1), q, 2);
        
        % Barycentric coordinate v check
        checkV = (v >= 0.0) & (u + v <= 1.0);
        if ~any(checkV), continue; end
        
        % Narrow down to triangles passing v-test
        idxFinal = idxStepU(checkV);
        f = f(checkV);
        q = q(checkV, :);
        
        t = f .* dot(e2(idxFinal, :), q, 2);
        
        % Check for positive intersection distance
        checkT = t > 1e-5;
        if ~any(checkT), continue; end
        
        % --- Identify Closest Intersection ---
        validT = t(checkT);
        validFinalIdx = idxFinal(checkT);
        [minT, localIdx] = min(validT);
        bestTriIdx = validFinalIdx(localIdx);
        
        % --- Record Intersection Result ---
        intersectPoints(r, :) = orig + minT * dir;
        
        % --- Reflection Calculation ---
        % Use the discrete face normal of the specific triangle hit
        nHit = faceNormals(bestTriIdx, :);
        
        % Ensure the normal vector points towards the source (handle back-face)
        if dot(dir, nHit) > 0
            nHit = -nHit;
        end
        
        % Law of Reflection: R = I - 2*(N·I)*N
        reflectDirs(r, :) = dir - 2 * dot(dir, nHit) * nHit;
    end
    
    % Final normalization of output directions
    reflectDirs = reflectDirs ./ vecnorm(reflectDirs, 2, 2);
end