function [r, dr] = ComputeEnrichedRadialBasis(point, neighs, sd_rad, options)
%% ComputeEnrichedRadialBasis
% Use: Computes an enriched radial basis and its gradient.
%
% Syntax: [r, dr] = ComputeEnrichedRadialBasis(point, neighs, sd_rad, options)
%
% Input:
%   point     - The coordinates of the evaluation point, format [1 x dim]: where dim the coordinates dimension
%   neighs    - The neighbor nodes indices of the point, format [k x 1], where k the number of neighbors
%   sd_rad    - The radius of the support domain, format: [1 x 1]
%   options   - The options for the generation of the meshfree approximation's weight function, format: struct
%
% Output:
%   r  - The value of enriched radial basis, format: [(k+m) x 1] where k the number of neighbors and m the number of
%        monomials in the linear monomial basis
%   dr - The gradient of the enriched radial basis, format: [(k+m) x dim] where k the number of neighbors and m the number of
%        monomials in the linear monomial basis
%
% Author: Konstantinos A. Mountris
% web: https://www.mountris.org
% mail: konstantinos.mountris@gmail.com
% license: see LICENSE.txt
%%      

% The number of neighbor nodes.
n = size(neighs,1);

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

% Initialize the enriched radial basis and its gradient.
r = zeros(n+m,1);
dr = zeros(n+m,dim);

% Iterate neighbor nodes.
for i = 1:n    
    % Compute weight function and gradient at current neighbor node.    
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

% Add monomial enrichment.
[p,dp] = MonomialBasis(point, m);
r(n+1:end) = p;
dr(n+1:end,:) = dp';
% for d = 2:dim+1
%     r(n+d) = point(d-1);
%     dr(n+d,d-1) = 1;
% end
       
end

