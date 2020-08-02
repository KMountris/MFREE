function [w, dw] = Multiquadric(point, sd_center, sd_rad, options)
%% Multiquadric
% Use: Computes the value of a multiquadric radial basis function function and its gradient at a given point.
%      The multiquadric radial basis function function is defined at the given support domain. 
%
% Syntax: [w, dw] = Multiquadric(point, sd_center, sd_rad, options)
%
% Input:
%   point     - The coordinates of the point, format [1 x dim]: where dim the coordinates dimension
%   sd_center - The coordinates of the support domain's center, format: [1 x dim] where dim the coordinates dimension
%   sd_rad    - The radius of the support domain, format: [1 x 1]
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

% Multiquadric RBF parameters.
beta = options.beta;
theta = options.theta;

% Multiquadric RBF Shape coefficient.
rc = theta*sd_rad;

% Compute the squared radius from the support domain center to the point.
r2 = sum((point-sd_center).^2);

% Compute the weight function and its gradient.
w = (r2 + rc*rc)^beta;

% 1/r * dw/dr
inv_r_dw_dr = 2*beta*(r2 + rc*rc)^(beta-1);      

% dw/dx = (1/r * dw/dr) * (dr/dx * r)
dw = inv_r_dw_dr .* (point-sd_center);

end

