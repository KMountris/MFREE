function [phi, dphi] = RpiShape(point, nodes, neighs, sd_rad, options)
%% RpiShape
% Use: Computes the meshfree Radial Point Interpolation (RPI) approximation basis function and its gradient.
%
% Syntax: [phi, dphi] = RpiShape(point, nodes, neighs, sd_rad, options)
%
% Input:
%   point     - The coordinates of the evaluation point, format [1 x dim]: where dim the coordinates dimension
%   nodes     - Node coordinates, format: [n x dim] where n the number of nodes and dim the coordinates dimension
%   neighs    - The neighbor nodes indices of the point, format [k x 1], where k the number of neighbors
%   sd_rad    - The radius of the support domain, format: [1 x 1]
%   options   - The options for the generation of the meshfree approximation's weight function, format: struct
%
% Output:
%   phi  - The value of the RPI basis function on the neighbor nodes, format: [k x 1] 
%          where k the number of neighbors
%   dphi - The gradient of the RPI basis function, format: [k x dim] 
%          where k the number of neighbors and dim the gradient dimension
%
% Author: Konstantinos A. Mountris
% web: https://www.mountris.org
% mail: konstantinos.mountris@gmail.com
% license: see LICENSE.txt
%%

% Number of the neighbor nodes.
neighs_num = size(neighs,1);

% The number of monimials in the monomial basis.
dim = size(point,2);
if strcmp(options.monomial,'linear')
    m = dim+1;
elseif strcmp(options.monomial,'quadratic')
    if dim == 1
        m = 3;
    elseif dim == 2
        m = 6;
    elseif dim == 3
        m = 10; 
    end
end

% Assemble the G matrix.
G = zeros(neighs_num+m,neighs_num+m);
for i = 1:neighs_num
    nn = neighs(i);
    
    rr = ComputeEnrichedRadialBasis(nodes(nn,:), nodes(neighs,:), sd_rad, options);
    
    G(:,i) = rr;
    p = MonomialBasis(nodes(nn,:), m);
    G(i,neighs_num+1:end) = p;
    G(neighs_num+1:end,i) = p;
end

% Compute RPI basis function and gradient.
[rr, drr] = ComputeEnrichedRadialBasis(point, nodes(neighs,:), sd_rad, options);

phi = G\rr;
phi = phi(1:neighs_num);

dphi = zeros(neighs_num+m, dim);
for d=1:dim
    dphi(:,d) = G\drr(:,d);
end
dphi = dphi(1:neighs_num,:);


end

