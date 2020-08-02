function [nodes, cells] = HexMesh(l, h)
%% HexMesh
% Use: Creates a regular hexahedral mesh in 3 dimensions.
%
% Syntax: [nodes, cells] = HexMesh(l, h)
%
% Input:
%   l - The edge size, format: [1 x 1]
%   h - The mesh spacing, format: [1 x 1]
%
% Output:
%   nodes - Node coordinates, format: [n x 3] where n the number of nodes
%   cells - Cells connectivity list, format [e x 8], where e the number of cells
%
% Author: Konstantinos A. Mountris
% web: https://www.mountris.org
% mail: konstantinos.mountris@gmail.com
% license: see LICENSE.txt
%%

% Create nodes.
[X,Y,Z] = meshgrid(0:h:l);
nodes = [X(:), Y(:), Z(:)];

% Nodes number distributed over the axis of the domain
nd_x = length(0:h:l);
nd_y = length(0:h:l);
nd_z = length(0:h:l);

% Find number of cells per axis of the domain and total
cell_x = nd_x-1;
cell_y = nd_y-1;
cell_z = nd_z-1;
cells_num = cell_x * cell_y * cell_z;

% Compute cells.
cells = zeros(cells_num,8);
for i = 1:cell_x
    for j = 1:cell_y
        for k = 1:cell_z
            id = cell_y*cell_z*(i-1) + cell_z*(j-1) + k;
            cells(id,:) = [(i-1)*nd_y*nd_z + (j-1)*nd_z + k, ...
                            (i-1)*nd_y*nd_z +  j*nd_z    + k, ...
                            (i-1)*nd_y*nd_z +  j*nd_z    + k+1, ...
                            (i-1)*nd_y*nd_z + (j-1)*nd_z + k+1, ...
                            i*nd_y*nd_z    + (j-1)*nd_z + k, ...
                            i*nd_y*nd_z    +  j*nd_z    + k, ...
                            i*nd_y*nd_z    +  j*nd_z    + k+1, ...
                            i*nd_y*nd_z    + (j-1)*nd_z + k+1]; 
        end
    end    
end

end

