function [nodes, cells] = QuadMesh(l, h)
%% QuadMesh
% Use: Creates a regular quadrilateral mesh in 2 dimensions.
%
% Syntax: [nodes, cells] = QuadMesh(l, h)
%
% Input:
%   l - The edge size, format: [1 x 1]
%   h - The mesh spacing, format: [1 x 1]
%
% Output:
%   nodes - Node coordinates, format: [n x 2] where n the number of nodes
%   cells - Cells connectivity list, format [e x 4], where e the number of cells
%
% Author: Konstantinos A. Mountris
% web: https://www.mountris.org
% mail: konstantinos.mountris@gmail.com
% license: see LICENSE.txt
%%

% Create nodes.
[X,Y] = meshgrid(0:h:l);
nodes = [X(:), Y(:)];

% Nodes number distributed over the length and the width of the domain
nd_row = length(0:h:l);
nd_col = length(0:h:l);

% Find number of cells per row per col and total
cell_row = nd_row-1;
cell_col = nd_col-1;
cells_num = cell_row * cell_col;

% Compute cells.
cells = zeros(cells_num,4);
for i = 1:cell_row
    for j = 1:cell_col    
        
        id = cell_col*(i-1) + j;
        cells(id,:) = [(i-1)*nd_col+j, i*nd_col+j,  i*nd_col+j+1, (i-1)*nd_col+j+1]; 
        
    end    
end

end

