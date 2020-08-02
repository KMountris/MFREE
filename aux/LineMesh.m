function [nodes, cells] = LineMesh(l, h)
%% LineMesh
% Use: Creates a 1D line mesh
%
% Syntax: [nodes, cells] = LineMesh(l, h)
%
% Input:
%   l - The edge size, format: [1 x 1]
%   h - The mesh spacing, format: [1 x 1]
%
% Output:
%   nodes - Node coordinates, format: [n x 1] where n the number of nodes
%   cells - Cells connectivity list, format [e x 2], where e the number of cells
%
% Author: Konstantinos A. Mountris
% web: https://www.mountris.org
% mail: konstantinos.mountris@gmail.com
% license: see LICENSE.txt
%%


% Create nodes.
nodes = (0:h(1):l(1))';

% Find number of cells
cells_num = length(nodes) -1;

% Compute cells.
cells = zeros(cells_num,2);
for id = 1:cells_num  
        cells(id,:) = [id, id+1]; 
end


end

