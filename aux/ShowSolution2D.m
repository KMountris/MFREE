function ShowSolution2D(nodes, cells, sol)
%% ShowSolution2D
% Use: Creates a patch plot to display a solution vector in 2 dimensions.
%
% Syntax: ShowSolution2D(nodes, cells, sol)
%
% Input:
%   nodes     - Node coordinates, format: [n x dim] where n the number of nodes and dim the coordinates dimension
%   cells     - Cells connectivity list, format: [e x k] where e the number of cells and k the number of nodes in the cell
%   sd_rad    - The radius of the support domain, format: [1 x 1]
%
% Author: Konstantinos A. Mountris
% web: https://www.mountris.org
% mail: konstantinos.mountris@gmail.com
% license: see LICENSE.txt
%%
nodes = [nodes(:,1), nodes(:,2), sol];    
patch('Vertices',nodes,'Faces',cells,'FaceColor','flat', 'FaceVertexCData', sol);


end