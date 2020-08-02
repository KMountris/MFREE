function [nodes, cells] = TriMesh(l, h)
%% TriMesh
% Use: Creates a 2D triangular mesh
%
% Syntax: [nodes, cells] = TriMesh(l, h)
%
% Input:
%   l - The edge size, format: [1 x 1]
%   h - The mesh spacing, format: [1 x 1]
%
% Output:
%   nodes - Node coordinates, format: [n x 2] where n the number of nodes
%   cells - Cells connectivity list, format [e x 3], where e the number of cells
%
% Author: Konstantinos A. Mountris
% web: https://www.mountris.org
% mail: konstantinos.mountris@gmail.com
% license: see LICENSE.txt
%%

if length(l) < 2
    L1 = l;
    L2 = l;
else
    L1 = l(1);
    L2 = l(2);
end

if length(h) < 2
    h = h;
else
    h = max(h(1),h(2));
end

% Create square geometry and mesh with 3 cm side.
gd = [ 3 ;   4; 
       0 ; L1;
      L1 ;  0;
       0 ;  0;
      L2 ; L2];
  
[geom] = decsg(gd);

[p,~,t] = initmesh(geom,'Hmax',h);

% Mesh data.
nodes = p(1:2,:)';
cells = t(1:3,:)';

end

