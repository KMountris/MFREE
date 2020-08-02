function [w, dw] = CubicSpline(point, sd_center, sd_rad, options)
%% CubicSpline
% Use: Computes the value of a cubic spline weight function and its gradient at a given point.
%      The cubic spline is defined at the given support domain. 
%
% Syntax: [w, dw] = CubicSpline(point, sd_center, sd_rad, options)
%
% Input:
%   point     - The coordinates of the point, format [1 x dim]: where dim the coordinates dimension
%   sd_center - The coordinates of the support domain's center, format: [1 x dim] where dim the coordinates dimension
%   sd_rad    - The radius of the support domain, format: [1 x 1]
%   options   - The options for the weight function generation, format: struct
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
r = sqrt(sum((point-sd_center).^2)) / (dil_sd_rad);

% Compute the weight function and its gradient.
w = 0.;
dw = zeros(1,dim);
if (r > 0.5 && r <= 1) 
    w = 4/3 - 4*r + 4*r*r - 4/3*r*r*r;
    % 1/r * dw/dr
    inv_r_dw_dr = -4/r + 8 - 4*r;
    
    % dw/dx = (1/r * dw/dr) * (dr/dx * r)
    dw = inv_r_dw_dr .* (point-sd_center) ./ (dil_sd_rad*dil_sd_rad);
elseif (r <= 0.5)
    w = 2/3 - 4*r*r + 4*r*r*r;
    
    % 1/r * dw/dr
    inv_r_dw_dr = -8 + 12*r;
    
    % dw/dx = (1/r * dw/dr) * (dr/dx * r)
    dw = inv_r_dw_dr .* (point-sd_center) ./ (dil_sd_rad*dil_sd_rad);
end


end
