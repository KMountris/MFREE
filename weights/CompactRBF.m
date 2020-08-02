function [w, dw] = CompactRBF(point, sd_center, sd_rad, options)
%% CompactRBF
% Use: Computes the value of a compact radial basis function and its gradient at a given point.
%      The compact radial basis function is defined at the given support domain. 
%
% Syntax: [w, dw] = CompactRBF(point, sd_center, sd_rad, options)
%
% Input:
%   point     - The coordinates of the point, format [1 x dim]: where dim the coordinates dimension
%   sd_center - The coordinates of the support domain's center, format: [1 x dim] where dim the coordinates dimension
%   sd_rad    - The radius of the support domain, format: [1 x 1]
%   options   - options
%
% Output:
%   w  - The value of the weight function, format: [1 x 1]
%   dw - The gradient of the weight function, format: [1 x dim] where dim the gradient dimension
%
% Author: Konstantinos A. Mountris
% web: https://www.mountris.org
% mail: konstantinos.mountris@gmail.com
% license: see LICENSE.txt
%%

% Coordinate dimension.
dim = size(point,2);

% Dilated support radius.
dil_sd_rad = options.dc*sd_rad;

% Compute the normalized radius from the support domain center to the point.
r = sqrt(sum((point-sd_center).^2)) / dil_sd_rad;

% Compute the weight function and its gradient.
w = 0.;
dw = zeros(1,dim);
if (r <= 1)
	w = (1-r)^6*(6 + 36*r + 82*r*r + 72*r*r*r + 30*r*r*r*r + 5*r*r*r*r*r);
   
    % 1/r * dw/dr
    inv_r_dw_dr = (11/dil_sd_rad)*(r-1)^5*(8 + 40*r + 48*r*r + 25*r*r*r + 5*r*r*r*r);
    
    % dw/dx = (1/r * dw/dr) * (dr/dx * r)
    dw = inv_r_dw_dr .* (point-sd_center) ./ dil_sd_rad;
end

end

