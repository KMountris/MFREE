function ShowSolution3D(points, ipoint, sol)
%% ShowSolution3D
% Use: Creates a slice plot to display a solution vector in 3 dimensions.
%
% Syntax: ShowSolution3D(points, ipoint, sol)
%
% Input:
%   points    - Points coordinates, format: [p x dim] where p the number of points and dim the coordinates dimension
%   ipoint    - Inspection point coordinates, format: [1 x dim] where dim the coordinates dimension
%   sol       - The solution vector, format: [p x 1] where p the number of points where the solution was calculated
%
% Author: Konstantinos A. Mountris
% web: https://www.mountris.org
% mail: konstantinos.mountris@gmail.com
% license: see LICENSE.txt
%%

n = int32(size(points,1)^(1/3));

eval_xmat = reshape(points(:,1),n,n,n);
eval_ymat = reshape(points(:,2),n,n,n);
eval_zmat = reshape(points(:,3),n,n,n);
sol_mat = reshape(sol,n,n,n);

slice(eval_xmat, eval_ymat, eval_zmat, sol_mat, ipoint(1), ipoint(2), ipoint(3));
  
end

