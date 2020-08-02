function AddMfreePaths
%% AddMfreePaths
% Use: Includes the file paths of MFREE in the MATLAB workspace. 
%      It should be called before using the MFREE functions. 
%
% Syntax: AddMfreePaths
%
% Author: Konstantinos A. Mountris
% web: https://www.mountris.org
% mail: konstantinos.mountris@gmail.com
% license: see LICENSE.txt
%%

% Add MFREE paths.
str = which('AddMfreePaths.m');
path = fileparts(str);
addpath(path);
addpath(strcat(path,'/aux'));
addpath(strcat(path,'/mfree-1.0.0'));
addpath(strcat(path,'/support'));
addpath(strcat(path,'/textprogressbar'));
addpath(strcat(path,'/weights'));

end
