function [p, dp] = MonomialBasis(point, m)
%% MonomialBasis
% Use: Computes the monomial basis for a given point.
% 
% Syntax: [p, dp] = MonomialBasis(point, m)
%
% Input:
%   point - The coordinates of the point, format [1 x dim]: where dim the coordinates dimension
%   m     - The number of the monomials in the monomial basis, format: [1 x 1]
%
% Output:
%   p  - The monomial basis for the given point, format: [1 x m]
%   dp - The gradient of the monomial basis, format: [dim x m] where dim the gradient dimension
%
% Author: Konstantinos A. Mountris
% web: https://www.mountris.org
% mail: konstantinos.mountris@gmail.com
% license: see LICENSE.txt
%%
dim = size(point,2);

x = point(1);
if dim > 1, y = point(2); end
if dim > 2, z = point(3); end
if dim > 3, error('Not supported dimension for monomial basis computation. Supported: 1 | 2 | 3 dimensions'); end

% Compute monomial basis for different order and spatial dimension
if m == dim+1
    if dim == 1
         p = [1 x];
        dp = [0 1];
    elseif dim == 2
         p = [1 x y];
        dp = [0 1 0;
              0 0 1];
    else
         p = [1 x y z];
        dp = [0 1 0 0;
              0 0 1 0;
              0 0 0 1];
    end
    
elseif m == 3 && dim == 1
     p = [1 x x*x];
    dp = [0 1 x];
elseif m == 6 && dim == 2
     p = [1 x y x*x y*y x*y];
    dp = [0 1 0 2*x  0  y;
          0 0 1  0  2*y x];
elseif m == 10 && dim == 3
     p = [1 x y z x*x y*y z*z x*y y*z x*z];
    dp = [0 1 0 0 2*x  0   0  y   0  z;
          0 0 1 0  0  2*y  0  x   z  0;
          0 0 0 1  0   0  2*z  0  y  x];
else
    error('Not supported monomial basis order');
end
    
end

