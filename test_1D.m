%% Test meshfree approximations for 1D domains
%
%  Inspection of the basis function and its first derivative
%  on a 1D cable using Moving Least Squares (MLS), Radial Point Interpolation (RPI), 
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
l = 10;     % length
h = 0.5;    % node spacing

% Geometry nodes.
[nodes, elem] = LineMesh(l,h);
nodes_num = length(nodes);

% Evaluation points.
eval = (0:h/10:l)';
eval_num = length(eval);

% Inspection point coordinates and index of nearest geometry node.
ipoint = l/2;
[~,inode] = min(dist(ipoint,nodes'));

% Determine support domain radius for each node of the geometry.
sd = SupportRadius(nodes, elem);

% Set neighbor nodes for each evaluation point. 
% We consider all the geometry nodes as neighbors.
neighs = (1:length(nodes))';

% Set options for MLS approximation.
mls_opt = MfreeOptions;
mls_opt.weight = 1;
mls_opt.dc = 1.8;

% Set options for RPI approximation.
rpi_opt = MfreeOptions;
rpi_opt.weight = 4;
rpi_opt.theta = 1.03;
rpi_opt.beta = 1.42;

% Set options for MKI approximation.
mki_opt = MfreeOptions;
mki_opt.weight = 3;
mki_opt.theta = 1.5;

% Set matrices to store basis functions and derivatives
MLS_PHI = zeros(eval_num, nodes_num);
MLS_PHIDX = zeros(eval_num, nodes_num);
RPI_PHI = zeros(eval_num, nodes_num);
RPI_PHIDX = zeros(eval_num, nodes_num);
MKI_PHI = zeros(eval_num, nodes_num);
MKI_PHIDX = zeros(eval_num, nodes_num);

% Compute basis function and derivative for each evaluation point.
textprogressbar('Calculating 1D basis function and gradient: ');
for i = 1:eval_num
    [mls_phi, mls_dphi] = MfreeShape(eval(i), nodes, neighs, sd, mls_opt, 'MLS');
    [rpi_phi, rpi_dphi] = MfreeShape(eval(i), nodes, neighs, sd, rpi_opt, 'RPI');
    [mki_phi, mki_dphi] = MfreeShape(eval(i), nodes, neighs, sd, mki_opt, 'MKI');
    
    % Collect basis function values.
    MLS_PHI(i,:) = mls_phi;
    RPI_PHI(i,:) = rpi_phi;
    MKI_PHI(i,:) = mki_phi;
    
    % Collect basis function derivative values.
    MLS_PHIDX(i,:) = mls_dphi;
    RPI_PHIDX(i,:) = rpi_dphi;
    MKI_PHIDX(i,:) = mki_dphi;
    
    j = (i/eval_num)*100;
    textprogressbar(j);
end
textprogressbar('done');

% Plot basis function.
figure; hold on;
plot(nodes,zeros(nodes_num,1),'o-.k','MarkerSize',8,'MarkerFaceColor','k');
plot(nodes(inode),0,'sk','MarkerSize',12,'MarkerFaceColor','m');
plot(eval,MLS_PHI(:,inode),'-b','LineWidth',3);
plot(eval,RPI_PHI(:,inode),'-.r','LineWidth',3);
plot(eval,MKI_PHI(:,inode),':g','LineWidth',3);
set(gca,'FontSize',16)
lgd = legend('points i','point I','MLS','RPI','MKI'); 
lgd.FontSize = 16;
% title('Basis function');
hold off;

% Plot X derivative of basis function.
figure; hold on;
plot(nodes,zeros(nodes_num,1),'o-.k','MarkerSize',8,'MarkerFaceColor','k');
plot(nodes(inode),0,'sk','MarkerSize',12,'MarkerFaceColor','m');
plot(eval, MLS_PHIDX(:,inode),'-b','LineWidth',3);
plot(eval, RPI_PHIDX(:,inode),'-.r','LineWidth',3);
plot(eval, MKI_PHIDX(:,inode),':g','LineWidth',3);
set(gca,'FontSize',16)
lgd = legend('points i','point I','MLS','RPI','MKI'); 
lgd.FontSize = 16;
% title('Basis function X derivative');
hold off;

