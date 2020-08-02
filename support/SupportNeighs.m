function [neighs] = SupportNeighs(points, nodes, sd_rad)
%% SupportNeighs
% Use: Computes the neighbor nodes for each evaluation point in a ND regular grid.
%
% Syntax: [neighs] = SupportNeighs(nodes, sd_radius)
%
% Input:
%   points - Points coordinates, format: [n x dim] where p the number of points and dim the coordinates dimension 
%   nodes  - Nodes coordinates, format: [n x dim] where n the number of nodes and dim the coordinates dimension
%   sd_rad - The support domain radius
%
% Output:
%   radius - The support domain radius for each node, format: [n x 1] where n the number of nodes
%
% Author: Konstantinos A. Mountris
% web: https://www.mountris.org
% mail: konstantinos.mountris@gmail.com
% license: see LICENSE.txt
%%

np = size(points,1);

%% Find neighbours of nodes
neighs = {[]};

% Iterate over the field nodes of the geometry.
for i = 1:np
    
    dist = sqrt(sum((nodes-points(i,:)).^2,2));
    neighs{i} = find(dist <= sd_rad);
end

end

