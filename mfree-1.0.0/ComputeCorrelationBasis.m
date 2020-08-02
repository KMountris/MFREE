function [r, dr] = ComputeCorrelationBasis(point, neighs, sd_rad, options)
%% ComputeCorrelationBasis
% Use: Computes the correlation basis and its gradient for the Moving Kriging Interpolation.
%
% Syntax: [r, dr] = ComputeCorrelationBasis(point, neighs, sd_rad, options)
%
% Input:
%   point     - The coordinates of the evaluation point, format [1 x dim]: where dim the coordinates dimension
%   neighs    - The neighbor nodes indices of the point, format [k x 1], where k the number of neighbors
%   sd_rad    - The radius of the support domain, format: [1 x 1]
%   options   - The options for the generation of the meshfree approximation's weight function, format: struct
%
% Output:
%   r  - The value of correlation basis, format: [k x 1] where k the number of neighbors
%   dr - The gradient of the correlation basis, format: [k x dim] where k the number of neighbors
%
% Author: Konstantinos A. Mountris
% web: https://www.mountris.org
% mail: konstantinos.mountris@gmail.com
% license: see LICENSE.txt
%%     

% The number of neighbor nodes.
neighs_num = size(neighs,1);

dim = size(point,2);

% Initialize the correlation basis and its gradient.
r = zeros(neighs_num,1);
dr = zeros(neighs_num,dim);

% Iterate neighbor nodes.
for i = 1:neighs_num

    % Compute weight function and its gradient for current neighbor.    
    if options.weight == 1
        [wi, dwi] = CubicSpline(point, neighs(i,:), sd_rad(i), options);
    elseif options.weight == 2
        [wi, dwi] = QuarticSpline(point, neighs(i,:), sd_rad(i), options);
    elseif options.weight == 3
        [wi, dwi] = Gaussian(point, neighs(i,:), sd_rad(i), options);
    elseif options.weight == 4
        [wi, dwi] = Multiquadric(point, neighs(i,:), sd_rad(i), options); 
    elseif options.weight == 5
        [wi, dwi] = ThinPlateSpline(point, neighs(i,:), options);
    elseif options.weight == 6
        [wi, dwi] = CompactRBF(point, neighs(i,:), sd_rad(i), options);
    else
        error('Unknown weight function type');
    end
    
    r(i) = wi;
    dr(i,:) = dwi;

end
       
end

