function [w, dw] = Gaussian(point, sd_center, sd_rad, options)
%% Gaussian
% Use: Computes the value of a Gaussian weight function and its gradient at a given point.
%      The Gaussian is defined at the given support domain. 
% 
% Syntax: [w, dw] = Gaussian(point, sd_center, sd_rad, options)
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

% Compute the radius from the support domain center to the point.
r = sqrt(sum((point-sd_center).^2));

% Compute shape parameter.
theta = options.theta;
a = theta / (sd_rad*sd_rad);

% Compute the weight function and its gradient.
% w = 0.;
% dw = zeros(1,dim);
% if (r <= 1) 

    w = exp(-a*r*r);
    % 1/r * dw/dr
    inv_r_dw_dr = -2*a*w;
      
    % dw/dx = (1/r * dw/dr) * (dr/dx * r)
    dw = inv_r_dw_dr .* (point-sd_center);
% end


end



