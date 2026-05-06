%% Get Ground Truth for a Dual-OAP System via Analytic Quadratic Solution
function [P2, P3, D3] = get_oap_system_ground_truth_matrix(P1, D1, V1, F1, V2, F2)
% GET_OAP_SYSTEM_GROUND_TRUTH_MATRIX Calculates analytic ray tracing for two OAPs
%
% USAGE:
%   [P2, P3, D3] = get_oap_system_ground_truth_matrix(P1, D1, V1, F1, V2, F2)
%
% INPUTS:
%   P1 - Nx3 incident ray positions
%   D1 - Nx3 incident ray directions
%   V1, F1 - Vertex and Focus of the first OAP
%   V2, F2 - Vertex and Focus of the second OAP
%
% OUTPUTS:
%   P2 - Nx3 intersection points on the first OAP
%   P3 - Nx3 intersection points on the second OAP
%   D3 - Nx3 final reflected ray directions
%
% DESCRIPTION:
%   This function provides the analytic "ground truth" solution by solving
%   the quadratic intersection equation for ideal parabolic surfaces. 
%   It processes all rays simultaneously using vectorized matrix operations.

    % 1. Ensure input alignment (Nx3)
    if size(P1, 2) ~= 3, P1 = P1'; end
    if size(D1, 2) ~= 3, D1 = D1'; end
    
    % 2. Pre-processing
    currentP = P1; 
    currentD = D1 ./ sqrt(sum(D1.^2, 2)); % Normalize incident directions
    
    % Store vertices and foci as 1x3 row vectors
    vList = {V1(:)', V2(:)'}; 
    fList = {F1(:)', F2(:)'};
    hitPoints = cell(1, 2);

    for i = 1:2
        % Current mirror parameters
        vCurr = vList{i}; fCurr = fList{i};
        fLength = norm(fCurr - vCurr);
        axisDir = (fCurr - vCurr) / fLength; % Optical axis (1x3)
        
        % 3. Prepare quadratic equation coefficients: At^2 + Bt + C = 0
        W = currentP - fCurr;         % Nx3 relative to focus
        V_mat = currentP - vCurr;     % Nx3 relative to vertex
        
        % Fast dot product via matrix multiplication/vectorization
        cosTheta = currentD * axisDir';                
        k_val = V_mat * axisDir' + fLength;            
        
        % Analytic coefficients (Nx1)
        A = 1 - cosTheta.^2;              
        B = 2 * (sum(currentD .* W, 2) - cosTheta .* k_val);  
        C = sum(W .* W, 2) - k_val.^2;      
        
        % 4. Robust root finding
        discriminant = B.^2 - 4.*A.*C;
        discriminant(discriminant < 0) = NaN; % Mark misses as NaN
        
        sqrtDelta = sqrt(discriminant);
        t_candidates = NaN(size(A, 1), 2);
        
        % Quadratic case
        isQuadratic = abs(A) > 1e-12;
        t_candidates(isQuadratic, 1) = (-B(isQuadratic) + sqrtDelta(isQuadratic)) ./ (2*A(isQuadratic));
        t_candidates(isQuadratic, 2) = (-B(isQuadratic) - sqrtDelta(isQuadratic)) ./ (2*A(isQuadratic));
        
        % Linear case (A approaches 0)
        t_candidates(~isQuadratic, 1) = -C(~isQuadratic) ./ B(~isQuadratic);
        
        % Select the smallest positive root (t > 1e-6)
        t_candidates(t_candidates <= 1e-6) = NaN;
        t = min(t_candidates, [], 2, 'omitnan'); 
        
        if any(isnan(t))
            warning('Mirror %d: %d rays missed the surface.', i, sum(isnan(t)));
        end
        
        % 5. Update intersection positions
        currentP = currentP + t .* currentD;
        hitPoints{i} = currentP;
        
        % 6. Calculate analytic reflection direction
        % Surface normal n = unit(unit(P-F) - axis)
        PF = currentP - fCurr;
        PF_norm = PF ./ sqrt(sum(PF.^2, 2));
        normalVec = PF_norm - axisDir; 
        normalVec = normalVec ./ sqrt(sum(normalVec.^2, 2));
        
        % Ensure normal faces the incident ray
        dotDN = sum(currentD .* normalVec, 2);
        flipMask = dotDN > 0;
        normalVec(flipMask, :) = -normalVec(flipMask, :);
        dotDN(flipMask) = -dotDN(flipMask);
        
        % Law of Reflection: D_out = D_in - 2*(D_in·n)*n
        currentD = currentD - 2 * dotDN .* normalVec;
        currentD = currentD ./ sqrt(sum(currentD.^2, 2)); 
    end

    % Assign outputs
    P2 = hitPoints{1}; 
    P3 = hitPoints{2}; 
    D3 = currentD;
end