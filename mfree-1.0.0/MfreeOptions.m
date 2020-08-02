function options = MfreeOptions()
%% MfreeOptions
% Use: Set default options for the mesh free approximation.
%
% Syntax: options = MfreeOptions()
%
% Output:
%   options - The options of the meshfree approximation, format: struct
%          options.weight - The weight function type, format: [1 x 1]
%            supported: 1 - cubic spline
%                       2 - quartic spline
%                       3 - gaussian function
%                       4 - multiquadric
%                       5 - thin plate spline
%                       6 - compact radial basis function
%        options.monomial - The type of monomial basis, format [string],
%            supported: linear | quadratic
%          options.dc     - The dilatation coefficient of the support domain, format: [1 x 1]
%          options.theta  - Shape parameter, format: [1 x 1]
%          options.beta   - Exponent parameter, format: [1 x 1]
%
% Author: Konstantinos A. Mountris
% web: https://www.mountris.org
% mail: konstantinos.mountris@gmail.com
% license: see LICENSE.txt
%%

% Default values.
options.weight = 1;     % cubic spline weight function
options.monomial = 'linear';
options.dc = 1;
options.theta = 1;
options.beta = 1; 

end

