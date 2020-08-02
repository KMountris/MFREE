%% Test meshfree approximations for 2D domains
%
%  Inspection of the basis function and its first derivative
%  on a 2D surface using Moving Least Squares (MLS), Radial Point Interpolation (RPI), 
%  and Moving Kriging Interpolation (MKI) approximations
%
%  Author: Konstantinos A. Mountris
%  web: https://www.mountris.org
%  mail: konstantinos.mountris@gmail.com
%  license: see LICENSE.txt
%%
clear; close all;

% Include paths.
AddMfreePaths;

% Geometry characteristics.
l = 3;              % length
h = 0.5;            % node spacing
irregular = 0;      % select irregular domain

if irregular == 1
    [nodes, nelem] = TriMesh(l, h); 
    [eval, eelem] = TriMesh(l, 0.2*h);
else
    [nodes, nelem] = QuadMesh(l, h);
    [eval, eelem] = QuadMesh(l, 0.2*h);
end
nodes_num = length(nodes);
eval_num = length(eval);

% Inspection point coordinates and index of nearest geometry node.
ipoint = [0.5*l, 0.5*l];
[~,inode] = min(dist(ipoint,nodes'));


% Determine support domain radius for each node of the geometry.
sd = SupportRadius(nodes, nelem);

dc = 2.5;

% Set neighbor nodes for each evaluation point. 
% We consider all the geometry nodes as neighbors.
% neighs = (1:length(nodes))';
% neighs_num = length(neighs);

% Set options for MLS approximation.
mls_opt = MfreeOptions;
mls_opt.weight = 1;
mls_opt.monomial = 'linear';
mls_opt.dc = dc;

% Set options for RPI approximation.
rpi_opt = MfreeOptions;
rpi_opt.weight = 4;
rpi_opt.theta = 1.03;
rpi_opt.monomial = 'linear';
rpi_opt.beta = 1.42;
rpi_opt.dc = dc;

% Set options for MKI approximation.
mki_opt = MfreeOptions;
mki_opt.weight = 3;
mki_opt.monomial = 'linear';
mki_opt.beta = 1.2;
mki_opt.dc = dc;

[neighs] = SupportNeighs(eval, nodes, dc.*sd);
neighs_num = size(nodes,1);

% Set matrices to store basis functions and derivatives
MLS_PHI = zeros(eval_num, neighs_num);
MLS_PHIDX = zeros(eval_num, neighs_num);
MLS_PHIDY = zeros(eval_num, neighs_num);
RPI_PHI = zeros(eval_num, neighs_num);
RPI_PHIDX = zeros(eval_num, neighs_num);
RPI_PHIDY = zeros(eval_num, neighs_num);
MKI_PHI = zeros(eval_num, neighs_num);
MKI_PHIDX = zeros(eval_num, neighs_num);
MKI_PHIDY = zeros(eval_num, neighs_num);

% Compute basis function and derivative for each evaluation point.
textprogressbar('Calculating 2D basis function and gradient: ');
for i = 1:eval_num
    [mls_phi, mls_dphi] = MfreeShape(eval(i,:), nodes, neighs{i}, sd, mls_opt, 'MLS');
    [rpi_phi, rpi_dphi] = MfreeShape(eval(i,:), nodes, neighs{i}, sd, rpi_opt, 'RPI');
    [mki_phi, mki_dphi] = MfreeShape(eval(i,:), nodes, neighs{i}, sd, mki_opt, 'MKI');
    
    % Collect basis function values.
    MLS_PHI(i,neighs{i}) = mls_phi;
    RPI_PHI(i,neighs{i}) = rpi_phi;
    MKI_PHI(i,neighs{i}) = mki_phi;
    
    % Collect basis function X derivative values.
    MLS_PHIDX(i,neighs{i}) = mls_dphi(:,1);
    RPI_PHIDX(i,neighs{i}) = rpi_dphi(:,1);
    MKI_PHIDX(i,neighs{i}) = mki_dphi(:,1);
    
    % Collect basis function Y derivative values.
    MLS_PHIDY(i,neighs{i}) = mls_dphi(:,2);
    RPI_PHIDY(i,neighs{i}) = rpi_dphi(:,2);
    MKI_PHIDY(i,neighs{i}) = mki_dphi(:,2);
    
    j = (i/eval_num)*100;
    textprogressbar(j);
end
textprogressbar('done');

% Plot basis function.
figure; hold on;
subplot(1,3,1);
ShowSolution2D(eval, eelem, MLS_PHI(:,inode))
title('MLS Basis function'); view([50,10]);
subplot(1,3,2);
ShowSolution2D(eval, eelem, RPI_PHI(:,inode))
title('RPI Basis function'); view([50,10]);
subplot(1,3,3);
ShowSolution2D(eval, eelem, MKI_PHI(:,inode))
title('MKI Basis function'); view([50,10]);

% Plot X derivative of basis function.
figure; hold on;
subplot(1,3,1);
ShowSolution2D(eval, eelem, MLS_PHIDX(:,inode))
title('MLS Basis function X derivative'); view([50,10]);
subplot(1,3,2);
ShowSolution2D(eval, eelem, RPI_PHIDX(:,inode))
title('RPI Basis function X derivative'); view([50,10]);
subplot(1,3,3);
ShowSolution2D(eval, eelem, MKI_PHIDX(:,inode))
title('MKI Basis function X derivative'); view([50,10]);

% Plot Y derivative of basis function.
figure; hold on;
subplot(1,3,1);
ShowSolution2D(eval, eelem, MLS_PHIDY(:,inode))
title('MLS Basis function Y derivative'); view([50,10]);
subplot(1,3,2);
ShowSolution2D(eval, eelem, RPI_PHIDY(:,inode))
title('RPI Basis function Y derivative'); view([50,10]);
subplot(1,3,3);
ShowSolution2D(eval, eelem, MKI_PHIDY(:,inode))
title('MKI Basis function Y derivative'); view([50,10]);



