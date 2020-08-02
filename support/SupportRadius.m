function radius = SupportRadius(nodes, elem)
%% SupportRadius
% Use: Computes the support radiues for each node in a ND grid.
%
% Syntax: radius = SupportRadius(nodes)
%
% Input:
%   nodes - Node coordinates, format: [n x dim] where n the number of nodes and dim the coordinates dimension
%   elem  - Connectivity of the grid's elements
%
% Output:
%   radius - The support domain radius for each node, format: [n x 1] where n the number of nodes
%
% Author: Konstantinos A. Mountris
% web: https://www.mountris.org
% mail: konstantinos.mountris@gmail.com
% license: see LICENSE.txt
%%

% Get the nodes number and coordinates dimension.
[n, dim] = size(nodes);

radius = zeros(n,1);
conn_nodes_num = zeros(n,1);

% Iterate over the cells of the grid.
for el = 1:size(elem,1)
    
    % Get the cell connected facets
    if dim == 1
        conn_facets = elem(el,:);
    elseif dim == 2
        conn_facets = ExtractEdges2D(elem(el,:));
    elseif dim == 3 && size(elem,2) == 4
        conn_facets = ExtractEdgesTet(elem(el,:));
    elseif dim == 3 && size(elem,2) == 8
        conn_facets = ExtractEdgesHex(elem(el,:));
    else 
        error('not supported element');
    end
    
    % Iterate over connected facets connectivity.
    for i = 1:size(conn_facets,1)
        % Iterate over nodes before the current one.
        dist = sqrt(sum((nodes(conn_facets(i,1),:) - nodes(conn_facets(i,2),:)).^2));

        % Add the distance value to the support of the influence node.
        radius(conn_facets(i,1)) = radius(conn_facets(i,1)) + dist;
        radius(conn_facets(i,2)) = radius(conn_facets(i,2)) + dist;
            
        % Accumulate number of connected nodes to the influence nodes of the element.
        conn_nodes_num(conn_facets(i,1)) = conn_nodes_num(conn_facets(i,1)) + 1;
        conn_nodes_num(conn_facets(i,2)) = conn_nodes_num(conn_facets(i,2)) + 1;
    end
    
%     figure; hold on;
%     patch('Vertices',nodes,'Faces',elem(el,:),'FaceColor','flat', 'FaceVertexCData', 1);
%     for i = 1:size(elem,2)
%         if size(nodes,2) == 2
%             text(nodes(elem(el,i),1),nodes(elem(el,i),2),num2str(i),'FontSize',15);
%         elseif size(nodes,2) == 3
%             text(nodes(elem(el,i),1),nodes(elem(el,i),2),nodes(elem(el,i),3),num2str(i),'FontSize',15);
%         end
%     end
%     hold off;
%     pause;
    
end

radius = radius./conn_nodes_num;


end




function edges = ExtractEdges2D(elem)

el_dim = size(elem,2);
edges = zeros(el_dim,2);

for i = 1:el_dim-1
    edges(i,:) = [elem(i), elem(i+1)];
end

edges(el_dim,:) = [elem(el_dim), elem(1)];

end


function edges = ExtractEdgesTet(elem)

el_dim = size(elem,2);
edges = zeros(6,2);

pos = 1;
for i = 1:3
    for j = i+1:4
        edges(pos,:) = [elem(i), elem(j)];
        pos = pos+1;
    end
end

end

function edges = ExtractEdgesHex(elem)

n1 = elem(1); n2 = elem(2); 
n3 = elem(3); n4 = elem(4);
n5 = elem(5); n6 = elem(6);
n7 = elem(7); n8 = elem(8);

edges = [n1 n2;
         n2 n3;
         n3 n4;
         n4 n1;
         n5 n6;
         n6 n7;
         n7 n8;
         n8 n5;
         n2 n6;
         n3 n7;
         n4 n8;
         n1 n5];

end
