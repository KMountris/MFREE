function [w, dw] = Weight(point, sd_center, sd_rad, options)
%% Weight
% Use: Computes the weight function and its gradient at a given point.
%      The weight function is defined at the given support domain. 
%
% Syntax: [w, dw] = Weight(point, sd_center, sd_rad, options)
%
% Input:
%   point     - The coordinates of the evaluation point, format [1 x dim]: where dim the coordinates dimension
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

if options.weight == 1
    [w, dw] = CubicSpline(point, sd_center, sd_rad, options);
elseif options.weight == 2
    [w, dw] = QuarticSpline(point, sd_center, sd_rad, options);
elseif options.weight == 3
    [w, dw] = Gaussian(point, sd_center, sd_rad, options);
elseif options.weight == 4
    [w, dw] = Multiquadric(point, sd_center, sd_rad, options); 
elseif options.weight == 5
    [w, dw] = ThinPlateSpline(point, sd_center, options); 
elseif options.weight == 6
    [w, dw] = CompactRBF(point, sd_center, sd_rad, options);
else
    error('Unknown weight function type');
end

end

