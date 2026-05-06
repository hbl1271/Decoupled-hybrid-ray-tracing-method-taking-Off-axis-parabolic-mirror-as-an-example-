%% Export Matrix Data to Zemax Grid Sag Format
function matrix2zemax(surfaceData, D, decenterX, decenterY, savePath)
% MATRIX2ZEMAX Converts a 2D surface matrix into a Zemax-readable .DAT file
%
% USAGE:
%   matrix2zemax(surfaceData, D, decenterX, decenterY, savePath)
%
% INPUTS:
%   surfaceData - 2D matrix (nx by ny) representing surface sag or error
%   D           - Full width/diameter of the data grid (used to calculate spacing)
%   decenterX   - Decenter of the grid in X (as required by Zemax header)
%   decenterY   - Decenter of the grid in Y (as required by Zemax header)
%   savePath    - Full path and filename for the output .DAT file
%
% DESCRIPTION:
%   This function writes the surface matrix into the Zemax 'Grid Sag' format.
%   The header follows the standard:
%   nx ny dx dy zero decenterX decenterY
%   Followed by the sag values in a single column.

    % 1. Calculate grid dimensions and spacing
    [nx, ny] = size(surfaceData);
    
    % dx and dy are calculated based on the full aperture width D
    dx = D / (nx - 1);
    dy = D / (ny - 1);
    
    % 2. Reshape surface data for file writing
    % Note: Zemax expects the data to be written in a specific order; 
    % the original logic uses a 1x1x3 reshape-like flow via fprintf directly on re_E.
    re_E = reshape(surfaceData, nx, ny, 1);
    
    % 3. File I/O
    fid = fopen(savePath, 'wt');
    if fid == -1
        error('Failed to open file for writing at: %s', savePath);
    end
    
    % Write the Zemax Grid Sag Header:
    % nx ny dx dy zero decenterX decenterY
    fprintf(fid, '%d %d %2.6f %2.6f %d %2.6f %2.6f\n', nx, ny, dx, dy, 0, decenterX, decenterY);
    
    % Write the sag data values
    % Note: fprintf with a single format specifier will iterate through the matrix
    fprintf(fid, '%2.10f\n', re_E);
    
    fclose(fid);
end