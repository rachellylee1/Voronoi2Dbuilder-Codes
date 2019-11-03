function [graindistribution] = Voronoi2Dbuilder(ngrains,bndwidth,mult,bins)
%% Setup Voronoi Grain Interior Points

%Change to micrometer scale
um = 1;
%GB width replicated on both sides
bndwidth = bndwidth/2;

width = 10;
height = 10;

X = [rand([1 ngrains])' rand([1 ngrains])']*10*mult;
% X = [(gallery('normaldata',[1 ngrains],0)*width)' (gallery('normaldata',[1 ngrains],1)*height)'];
X = unique(X, 'rows'); %No duplicate points

outsidex = 100*width*[-10*ones(1, 18), -10:10, 10*ones(1, 18), 10:-1:-10]';
outsidey = 100*height*[-9:10, 10*ones(1, 18), 10:-1:-10, -10*ones(1, 19)]';

[VX, VY] = voronoi([X(:,1); outsidex], [X(:,2); outsidey]);

dt = delaunayTriangulation([X; outsidex, outsidey]);
[V,R] = voronoiDiagram(dt);

%% Find Voronoi Boundary Points
% Define Voronoi Boundary
closedist = abs(X(1,1)-X(1,2))/15;
closeoffx_l = min(X(:,1))-closedist; closeoffx_u = max(X(:,1))+closedist;
closeoffy_l = min(X(:,2))-closedist; closeoffy_u = max(X(:,2))+closedist;

%Add something so that bndline doesn't hit a vertex(but highly unlikely)
bndlinex = [closeoffx_l closeoffx_l closeoffx_u closeoffx_u closeoffx_l]; 
bndliney = [closeoffy_l closeoffy_u closeoffy_u closeoffy_l closeoffy_l];

Xold = X;
Rold = R;
Vold = V;

%% Create new vertices for all edges intersecting the bounding box
X = Xold;
R = Rold;
V = Vold;
% Grains that are updated.
updated = [];
% Grains that should be deleted.
delgrains = [];
% Seeds that should be deleted.
delseeds = [];
% Iterate over every grain.
for igrain = 1:size(R, 1)
    rnew = [];
    nedges = length(R{igrain});
    deletecell = true;
    % Iterate over every grain edge.
    for iedge = 1:nedges
        % Find the edge vertex indices and coordinates.
        ip1 = R{igrain}(iedge);
        p1 = V(ip1, :);
        if isinf(p1(1))
            continue;
        end
        ip2 = R{igrain}(1);
        if iedge ~= nedges
            ip2 = R{igrain}(iedge+1);
        end
        p2 = V(ip2, :);
        if isinf(p2(1))
            continue;
        end
        % Find which vertices are inside the bounding box.
        p1in = inpolygon(p1(1), p1(2), bndlinex, bndliney);
        p2in = inpolygon(p2(1), p2(2), bndlinex, bndliney);
        % If p1 is outside, then the old vertex is equal to iedge.
        if ~p1in
            iout = ip1;
            pout = p1;
        end
        if ~p2in
            % If p1 and p2 are outside, then the whole edge is outside, so
            % skip it.
            if ~p1in
                % Add only p1.
                rnew(end+1) = ip1;
                continue;
            end
            % If only p2 is outside, then the old vertex is equal to the
            % next edge after iedge.
            iout = ip2;
            pout = p2;
        end
        % No longer want to delete this grain, because at least one vertex
        % is inside the bounding box.
        deletecell = false;
        % If both p1 and p2 are inside, then skip this edge.
        if p1in && p2in
            % Add only p1.
            rnew(end+1) = ip1;
            continue;
        end
        % Find points of intersection between this edge and the bounding
        % box.
        [xi, yi] = polyxpoly([p1(1), p2(1)], [p1(2), p2(2)], bndlinex, bndliney, 'unique');
        % If no intersections, then something went wrong.
        if isempty(xi) || isempty(yi)
            fprintf('Error: No intersection\n');
            continue;
        end
        % If multiple intersections, choose the one closest to the point
        % outside the bounding box.
        if length(xi) > 1
            mindist = inf;
            xintersect = 0;
            for xx = 1:length(xi)
                testdist = sum((pout - [xi(xx), yi(xx)]).^2);
                if testdist < mindist
                    mindist = testdist;
                    xintersect = xx;
                end
            end
            pnew = [xi(xintersect), yi(xintersect)];
        else
            pnew = [xi, yi];
        end
        % Create a new vertex at the intersection point.
        ipnew = size(V, 1) + 1;
        V(ipnew, 1:2) = pnew;
        % Update this cell's vertices to include the new vertex.
        rnew(end+1) = ip1;
        rnew(end+1) = ipnew;
    end
    % Add to deleted grains
    if deletecell == true
        if igrain <= size(X, 1)
            delseeds(end+1) = igrain;
        end
        delgrains(end+1) = igrain;
    end
    % Update R-indices
    R{igrain} = rnew;
end
% Delete unused grains.
R(delgrains) = [];
X(delseeds, :) = [];
%% Make a list of duplicate verties.
dupes = [];
for ivert = 1:size(V, 1)
    % If this current vertex is a duplicate of a previous one, skip it.
    if any(dupes == ivert)
        continue;
    end
    % Test all cells to see if they have a point with identical
    % coordiantes. If so, add that point to the list of duplicates.
    for igrain = 1:length(R)
        for ipoint = 1:length(R{igrain})
            itest = R{igrain}(ipoint);
            if itest == ivert
                continue;
            end
            if V(itest, :) == V(ivert, :)
                R{igrain}(ipoint) = ivert;
                dupes(end+1) = itest;
            end
        end
    end
end
dupes = sort(unique(dupes), 'descend');
% Decrement all indices before deleting the duplicates.
for idupe = dupes
    for igrain = 1:length(R)
        for ipoint = 1:length(R{igrain})
            itest = R{igrain}(ipoint);
            if itest > idupe
                R{igrain}(ipoint) = itest - 1;
            end
        end
    end
end
% Finally, delete the duplicates.
V(dupes, :) = [];
%% Make a list of all points outside the bounding box.
outs = [];
tol = 1e-12;
for ivert = 1:size(V, 1)
    in = true;
    % Allow some tolerance, or else vertices on the bounding box might be
    % considered outside.
    if min(bndlinex) - V(ivert, 1) > tol || V(ivert, 1) - max(bndlinex) > tol
        in = false;
    end
    if min(bndliney) - V(ivert, 2) > tol || V(ivert, 2) - max(bndliney) > tol
        in = false;
    end
    if ~in
        outs(end+1) = ivert;
    end
end
outs = sort(unique(outs), 'descend');
% Decrement all indices before deleting the outside points.
for iout = outs
    for igrain = 1:length(R)
        removed = [];
        for ipoint = 1:length(R{igrain})
            itest = R{igrain}(ipoint);
            if itest == iout
                removed(end+1) = ipoint;
            elseif itest > iout
                R{igrain}(ipoint) = itest - 1;
            end
        end
        R{igrain}(removed) = [];
    end
end
% Finally, delete all outside points.
V(outs, :) = [];
%% Add corner vertices
corners = [ closeoffx_l, closeoffy_l;
            closeoffx_l, closeoffy_u;
            closeoffx_u, closeoffy_l;
            closeoffx_u, closeoffy_u];
for icorner = 1:size(corners, 1)
    if ~any(ismember(V, corners(icorner, :), 'rows'))
        % Search for cell with the closest two vertices
        mindist = inf;
        cornercell = [];
        for igrain = 1:length(R)
            % Find closest seed
            d = (X(igrain, 1) - corners(icorner, 1))^2 + ...
                (X(igrain, 2) - corners(icorner, 2))^2;
            if d < mindist
                mindist = d;
                cornercell = igrain;
            end
        end
        if isempty(cornercell)
            fprintf('Error: No corner cell found for: %.f\n', icorner);
            continue;
        end
        % Find the two vertices that do not intersect the cell when a line
        % is drawn to the corner vertex.
        nvert = length(R{cornercell});
        v1 = [];
        v2 = [];
        cellx = [V(R{cornercell}, 1); V(R{cornercell}(1), 1)];
        celly = [V(R{cornercell}, 2); V(R{cornercell}(1), 2)];
        for ivert = 1:nvert
            linex = [V(R{cornercell}(ivert), 1), corners(icorner, 1)];
            liney = [V(R{cornercell}(ivert), 2), corners(icorner, 2)];
            [xi, yi] = polyxpoly(linex, liney, cellx, celly, 'unique');
            if length(xi) == 1 && length(yi) == 1
                if isempty(v1)
                    v1 = ivert;
                elseif isempty(v2)
                    v2 = ivert;
                else
                    figure();
                    plot(V(R{cornercell}, 1), V(R{cornercell}, 2), '.b');
                    fprintf('Error: 3 or more vertices near corner %.f\n', icorner);
                end
            end
        end
        % Create the corner vertex
        V(end+1, 1:2) = corners(icorner, :);
        inew = size(V, 1);
        % Place the corner between the two vertices
        if v2 - v1 == 1
            R{cornercell} = [R{cornercell}(1:v1), inew, R{cornercell}(v2:end)];
        else
            R{cornercell} = [R{cornercell}, inew];
        end
    end
end

%% Plot Original Voronoi
% figure();
% hold on;

daspect([1 1 1]);
for igrain = 1:length(R)
    nedges = length(R{igrain});
    for iedge = 1:nedges
        p1 = V(R{igrain}(iedge), :);
        if isinf(p1(1))
            continue;
        end
        p2 = V(R{igrain}(1), :);
        if iedge ~= nedges
            p2 = V(R{igrain}(iedge+1), :);
        end
        if isinf(p2(1))
            continue;
        end
%         plot([p1(1), p2(1)], [p1(2), p2(2)], '-x');
    end
end
% title('Original Voronoi Plot')

%%
% R = cell array of coordinates of every seed
% Each cell of R holds an array of indices of the vertex coordinates
% V = coordinates of each vertex
% X = coordinates of each seed
Rnew = cell(0);
Vnew = [];
for icell = 1:size(R, 1)
    pgon = R{icell};
    Rnew{icell, 1} = [];
    for ip = 1:size(pgon, 2)
        ipnext = ip + 1;
        ipprev = ip - 1;
        if ip == size(pgon, 2)
            ipnext = 1;
        elseif ip == 1
            ipprev = size(pgon, 2);
        end
        ivert = pgon(ip);
        inext = pgon(ipnext);
        iprev = pgon(ipprev);
        % Math
        % Edge next
        dx = V(inext, 1) - V(ivert, 1);
        dy = V(inext, 2) - V(ivert, 2);
        mult = bndwidth / sqrt(dx^2 + dy^2);
        vna = [V(ivert, 1) - mult*dy, V(ivert, 2) + mult*dx];
        vnb = [V(inext, 1) - mult*dy, V(inext, 2) + mult*dx];
        
        % Edge prev
        dx = V(ivert, 1) - V(iprev, 1);
        dy = V(ivert, 2) - V(iprev, 2);
        mult = bndwidth / sqrt(dx^2 + dy^2);
        vpa = [V(iprev, 1) - mult*dy, V(iprev, 2) + mult*dx];
        vpb = [V(ivert, 1) - mult*dy, V(ivert, 2) + mult*dx];
        
        % Intersection
        p = vna;    r = (vnb - vna);
        q = vpa;    s = (vpb - vpa);
        t = cross2D(q - p, s) / cross2D(r, s);
        u = cross2D(q - p, r) / cross2D(r, s);
        indexV = size(Vnew, 1) + 1;
        if cross2D(r, s) ~= 0 && t >= 0 && t <= 1 && u >= 0 && u <= 1
            pintersect = p + t*r;
            Vnew(indexV, 1:2) = pintersect;
        else
            fprintf('Failure: t = %f, u = %f\n', t, u);
            Vnew(indexV, 1:2) = vpb;
        end
        Rnew{icell, 1}(ip) = indexV;
    end
    % Plot
%     figure();
    pgnew = Rnew{icell, 1};
    plot(V([pgon, pgon(1)], 1), V([pgon, pgon(1)], 2), 'x');
    plot(Vnew([pgnew, pgnew(1)], 1), Vnew([pgnew, pgnew(1)], 2), '-');
    drawnow;
    hold on;
    title('Grain & Grain Boundary Model')
    polyareas(icell) = polyarea(Vnew([pgnew,pgnew(1)],1), Vnew([pgnew,pgnew(1)],2));
    
%         Create text file for each Grain Interior Vertex points
    fileID = fopen(sprintf('G:\\Vertex\\Vertex%d.txt',icell),'w');
    fprintf(fileID,'%f %f\r\n',Vnew([pgnew pgnew(1)], :)');    
end

figure;
graindistribution = histogram(polyareas,bins);
title('Distribution of Grain Areas')
xlabel('Area (um^2)')
ylabel('Number of Grains')

%     polyareas(a) = polyarea(Vgrainb(:,1),Vgrainb(:,2));
%     % Plot New Voronoi w/ Grain Boundary
%     plot(Vgrainb(:,1),Vgrainb(:,2));
%     title('Grain & Grain Boundary Model')
%     drawnow;
%     % Create text file for each Grain Interior Vertex points
% %     fileID = fopen(sprintf('G:\\Vertex\\Vertex%d.txt',a),'w');
% %     fprintf(fileID,'%f %f\r\n',Vgrainb');
% end
% figure;
% graindistribution = histogram(polyareas,10);

% Create text file for Boundary around Voronoi
% fileID = fopen(sprintf('G:\\Vertex\\Boundary.txt'),'w');
% fprintf(fileID,'%f %f\r\n',[bndlinex;bndliney]);

