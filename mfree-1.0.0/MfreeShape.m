function [phi, dphi] = MfreeShape(point, nodes, neighs, sd_rad, options, type)
%% MfreeShape
% Use: Computes the meshfree approximation basis function and its gradient.
%
% Syntax: [phi, dphi] = MfreeShape(point, nodes, neighs, sd_rad, options, type)
%
% Input:
%   point     - The coordinates of the evaluation point, format [1 x dim]: where dim the coordinates dimension
%   nodes     - Node coordinates, format: [n x dim] where n the number of nodes and dim the coordinates dimension
%   neighs    - The neighbor nodes indices of the point, format [k x 1], where k the number of neighbors
%   sd_rad    - The radius of the support domain, format: [1 x 1]
%   options   - The options for the generation of the meshfree approximation's weight function, format: struct
%
% The type of the meshfree approximation depends on the value of the variable type
% Available types: MfreeType.MLS | MfreeType.RPI | MfreeType.MKI
%
% Output:
%   phi  - The value of the meshfree basis function on the neighbor nodes, format: [k x 1] 
%          where k the number of neighbors
%   dphi - The gradient of the meshfree basis function, format: [k x dim] 
%          where k the number of neighbors and dim the gradient dimension
%
% Author: Konstantinos A. Mountris
% web: https://www.mountris.org
% mail: konstantinos.mountris@gmail.com
% license: see LICENSE.txt
%%

if strcmp(type, 'MLS')
    [phi, dphi] = MlsShape(point, nodes, neighs, sd_rad, options);
elseif strcmp(type, 'RPI')
    [phi, dphi] = RpiShape(point, nodes, neighs, sd_rad, options);
elseif strcmp(type, 'MKI')
    [phi, dphi] = MkiShape(point, nodes, neighs, sd_rad, options); 
else
    error('Unknown meshfree approximation type');
end


end

