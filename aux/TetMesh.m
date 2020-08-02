function [nodes, cells] = TetMesh(l, h)
%% TetMesh
% Use: Creates a 3D tetrahedral mesh
%
% Syntax: [nodes, cells] = TetMesh(l, h)
%
% Input:
%   l - The edge size, format: [1 x 1]
%   h - The mesh spacing, format: [1 x 1]
%
% Output:
%   nodes - Node coordinates, format: [n x 3] where n the number of nodes
%   cells - Cells connectivity list, format [e x 4], where e the number of cells
%
% Author: Konstantinos A. Mountris
% web: https://www.mountris.org
% mail: konstantinos.mountris@gmail.com
% license: see LICENSE.txt
%%

L1 = l(1);
L2 = l(1);
L3 = l(1);

H1 = h(1);
H2 = h(1);
H3 = h(1);

% Create nodes.
[X,Y,Z] = meshgrid(0:H1:L1,0:H2:L2,0:H3:L3);
nodes = [X(:), Y(:), Z(:)];

% Nodes number distributed over the axis of the domain
nx = length(0:H1:L1);
ny = length(0:H2:L2);
nz = length(0:H3:L3);

% meshgrid flips x and y ordering
idx = reshape(1:prod([nx,ny,nz]),[nx,ny,nz]);
v1 = idx(1:end-1,1:end-1,1:end-1); v1=v1(:);
v2 = idx(1:end-1,2:end,1:end-1);   v2=v2(:);
v3 = idx(2:end,1:end-1,1:end-1);   v3=v3(:);
v4 = idx(2:end,2:end,1:end-1);     v4=v4(:);
v5 = idx(1:end-1,1:end-1,2:end);   v5=v5(:);
v6 = idx(1:end-1,2:end,2:end);     v6=v6(:);
v7 = idx(2:end,1:end-1,2:end);     v7=v7(:);
v8 = idx(2:end,2:end,2:end);       v8=v8(:);

cells = [v1  v3  v8  v7;
         v1  v8  v5  v7;
         v1  v3  v4  v8;
         v1  v4  v2  v8;
         v1  v6  v5  v8;
         v1  v2  v6  v8];

end

