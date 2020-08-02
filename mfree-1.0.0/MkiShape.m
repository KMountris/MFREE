function [phi, dphi] = MkiShape(point, nodes, neighs, sd_rad, options)
%% MkiShape
% Use: Computes the meshfree Moving Kriging Interpolation (MKI) approximation basis function and its gradient.
%
% Syntax: [phi, dphi] = RpiMkiShapeShape(point, nodes, neighs, sd_rad, options)
%
% Input:
%   point     - The coordinates of the evaluation point, format [1 x dim]: where dim the coordinates dimension
%   nodes     - Node coordinates, format: [n x dim] where n the number of nodes and dim the coordinates dimension
%   neighs    - The neighbor nodes indices of the point, format [k x 1], where k the number of neighbors
%   sd_rad    - The radius of the support domain, format: [k x 1], where k the number of neighbors
%   options   - The options for the generation of the meshfree approximation's weight function, format: struct
%
% Output:
%   phi  - The value of the MKI basis function on the neighbor nodes, format: [k x 1] 
%          where k the number of neighbors
%   dphi - The gradient of the MKI basis function, format: [k x dim] 
%          where k the number of neighbors and dim the gradient dimension
%
% Author: Konstantinos A. Mountris
% web: https://www.mountris.org
% mail: konstantinos.mountris@gmail.com
% license: see LICENSE.txt
%%

% Get the number of the support nodes.
neighs_num = size(neighs,1);

% Number of monomials in monomial basis.
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

% Initialize polynomial moment matrix.
P = zeros(neighs_num, m);

% ** Compute the R moment matrix
R = zeros(neighs_num, neighs_num);
for i = 1:neighs_num
    nn = neighs(i);
    
    gamma = ComputeCorrelationBasis(nodes(nn,:), nodes(neighs,:), sd_rad, options);
    
    R(i,:) = gamma;
    P(i,:) = MonomialBasis(nodes(nn,:), m);
end

% Compute inverse R and transpose P matrices.
I = eye(neighs_num);
Rinv = R\I;
Pt = P';

% Compute A and B matrices.
A = (Pt*Rinv*P)\(Pt*Rinv);
B = Rinv*(I-P*A);

% ***** Compute shape function and derivatives
[p, dp] = MonomialBasis(point,m);

[w, dw] = ComputeCorrelationBasis(point, nodes(neighs,:), sd_rad, options);

% MKI basis function
phi = p*A + w'*B;
phi = phi';

% MKI derivative
dphi = dp*A + dw'*B;
dphi = dphi';

end

