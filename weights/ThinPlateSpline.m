function [w, dw] = ThinPlateSpline(point, sd_center, options)
%% ThinPlate
% Use: Computes the value of a thin plate spline function and its gradient at a given point.
%      The thin plate spline function is defined at the given support domain. 
%
% Syntax: [w, dw] = ThinPlate(point, sd_center, options)
%
% Input:
%   point     - The coordinates of the point, format [1 x dim]: where dim the coordinates dimension
%   sd_center - The coordinates of the support domain's center, format: [1 x dim] where dim the coordinates dimension
%   options   - The options for the multiquadric radial basis function formation.
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

% Thin plate spline parameter.
beta = options.beta;

% Compute the radius from the support domain center to the point.
r = sqrt(sum((point-sd_center).^2));

% Compute the weight function and its gradient.
w = (r*r)^(0.5*beta);

% 1/r * dw/dr
inv_r_dw_dr = beta*(r*r)^(0.5*beta-1);

% dw/dx = (1/r * dw/dr) * (dr/dx * r)
dw = inv_r_dw_dr .* (point-sd_center);


end

