function [phi, dphi] = MlsShape(point, nodes, neigh_ids, sd_rad, options)
%% MlsShape
% Use: Computes the meshfree Moving Least Squares (MLS) approximation basis function and its gradient.
%
% Syntax: [phi, dphi] = MlsShape(point, nodes, neighs, sd_rad, options)
%
% Input:
%   point     - The coordinates of the evaluation point, format [1 x dim]: where dim the coordinates dimension
%   nodes     - Node coordinates, format: [n x dim] where n the number of nodes and dim the coordinates dimension
%   neighs    - The neighbor nodes indices of the point, format [k x 1], where k the number of neighbors
%   sd_rad    - The radius of the support domain, format: [1 x 1]
%   options   - The options for the generation of the meshfree approximation's weight function, format: struct
%
% Output:
%   phi  - The value of the MLS basis function on the neighbor nodes, format: [k x 1] 
%          where k the number of neighbors
%   dphi - The gradient of the MLS basis function, format: [k x dim] 
%          where k the number of neighbors and dim the gradient dimension
%
% Author: Konstantinos A. Mountris
% web: https://www.mountris.org
% mail: konstantinos.mountris@gmail.com
% license: see LICENSE.txt
%%

% Number of neighbor nodes.
n = size(neigh_ids,1);

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

% Initialize matrices
A = zeros(m,m);  dA = zeros(m,m,dim);   
B = zeros(m,n);  dB = zeros(m,n,dim);

% Initialize weight function
w  = zeros(n,1); 
dw = zeros(n,dim);

% Iterate neighbor nodes.
for i = 1:n
    
    neigh = nodes(neigh_ids(i),:);
    
    % Compute weight function value for current neigh.   
    if options.weight == 1
        [wi, dwi] = CubicSpline(point, neigh, sd_rad(i), options);
    elseif options.weight == 2
        [wi, dwi] = QuarticSpline(point, neigh, sd_rad(i), options);
    elseif options.weight == 3
        [wi, dwi] = Gaussian(point, neigh, sd_rad(i), options);
    elseif options.weight == 4
        [wi, dwi] = Multiquadric(point, neigh, sd_rad(i), options); 
    elseif options.weight == 5
        [wi, dwi] = ThinPlateSpline(point, neigh, options); 
    elseif options.weight == 6
        [wi, dwi] = CompactRBF(point, neigh, sd_rad(i), options);
    else
        error('Unknown weight function type');
    end
    
    [pxi, ~] = MonomialBasis(neigh, m);
    
    A = A + wi*(pxi'*pxi);
    B(:,i) = wi*pxi';
    
    for d = 1:dim
        dA(:,:,d) = dA(:,:,d) + dwi(d)*(pxi'*pxi);
        dB(:,i,d) = dwi(d)*pxi';
    end
    
end

% Correction for modified MLS with quadratic polynomial basis
if (m > dim+1)
    cf = 1e-7;
    for d = dim+2:m
        A(d,d) = A(d,d) + cf;
    end
end

% Compute MLS coefficients.
Ctrans = A\B;
C = Ctrans';

% Compute shape function and derivative.
[p, dp] = MonomialBasis(point, m);
phi = C*p';

dphi = zeros(n,dim);
for d = 1:dim
    B1 = dB(:,:,d) - dA(:,:,d)*Ctrans;
    
    coeffT1 = A\B1;
    coeff1 = coeffT1';
    
    dphi(:,d) = C*dp(d,:)' + coeff1*p';
    
end
    


end

