classdef Body
    properties
        numNodalPoints; NodalPoints;
        numTriangles; Triangles;
        numRectangles; Rectangles;
        Thickness;
        BoundaryEdges;
        numSubRegions; SubRegions;
        numContours; Contours;
        Density; lambda; mu; lambda_vis; mu_vis;
        strain_potential_energy;
        gravitational_potential_energy;
        J_lambda; J_mu;
        Stiffness_Matrix;
        Damping_Matrix;
        Inertia_Matrix;
        Gravitational_Vector;
        coef_order_one; coef_order_two; c2i; c2j; coef_order_three; c3i; c3j; c3k; % coefficient matrices for Green strain based forces
    end
    methods
        function obj = Body(npoints, points, ntris, tris, h)
            obj.numNodalPoints = npoints;
            for k=1:npoints
                pt(k) = NodalPoint(points(:,k));
            end
            obj.NodalPoints = pt;
            %
            if ntris > 0
                obj.numTriangles = ntris;
                for p=1:ntris
                    i = tris(p,1); j = tris(p,2); k = tris(p,3);
                    tr(p) = Triangle(i, j, k, points(:,i), points(:,j), points(:,k), h);
                end
                obj.Triangles = tr;
            end
            %
            obj.Thickness = h;
        
            obj.J_lambda = zeros(2*obj.numNodalPoints, 2*obj.numNodalPoints);
            obj.J_mu     = zeros(2*obj.numNodalPoints, 2*obj.numNodalPoints);
            for p=1:obj.numTriangles
                tri = obj.Triangles(p);
                vs = tri.Vertices;
                i = vs(1); j = vs(2); k = vs(3);
                loc = [ 2*i-1, 2*i, 2*j-1, 2*j, 2*k-1, 2*k ];
                obj.J_lambda(loc,loc) = obj.J_lambda(loc,loc) + tri.Partial_J_lambda;
                obj.J_mu(loc,loc)     = obj.J_mu(loc,loc)     + tri.Partial_J_mu;
            end

            edges = [];
            for p=1:obj.numTriangles
                tri = obj.Triangles(p);
                vs = tri.Vertices;
                i = vs(1); j = vs(2); k = vs(3);
                edges = [ edges; i, j; j, k; k, i ];
            end
            edges_sorted = sort(edges,2);
            n = size(edges_sorted, 1);
            onlyone = 1:n;
            for p=1:n-1
                [lia, loc] = ismember(edges_sorted(p,:), edges_sorted(p+1:n,:), 'rows');
                if lia
                    onlyone(p) = 0;
                    onlyone(p+loc) = 0;
                end
            end
            onlyone = setdiff(onlyone, 0);
            obj.BoundaryEdges = edges(onlyone, :);
            
            [ ncontours, contours, order ] = edges_to_contours (obj.BoundaryEdges);
            obj.numContours = ncontours;
            for p=1:ncontours
                [contour_p, area_p] = extract_contour(p, contours, order, obj.BoundaryEdges, points );
                ct(p) = Contour(size(contour_p,1), contour_p, area_p);
            end
            obj.Contours = ct;
        end

        function obj = rectangle_elements(obj, nrects, rects)
            h = obj.Thickness;
            %
            if nrects > 0
                obj.numRectangles = nrects;
                for p=1:nrects
                    i = rects(p,1); j = rects(p,2); k = rects(p,3); l = rects(p,4);
                    xi = obj.NodalPoints(i).Coordinates;
                    xj = obj.NodalPoints(j).Coordinates;
                    xl = obj.NodalPoints(l).Coordinates;
                    rc(p) = Rectangle(i, j, k, l, xj(1)-xi(1), xl(2)-xi(2), h);
                end
                obj.Rectangles = rc;
            end
            for p=1:obj.numRectangles
                rct = obj.Rectangles(p);
                vs = rct.Vertices;
                i = vs(1); j = vs(2); k = vs(3); l = vs(4);
                loc = [ 2*i-1, 2*i, 2*j-1, 2*j, 2*k-1, 2*k, 2*l-1, 2*l ];
                obj.J_lambda(loc,loc) = obj.J_lambda(loc,loc) + rct.Partial_J_lambda;
                obj.J_mu(loc,loc)     = obj.J_mu(loc,loc)     + rct.Partial_J_mu;
            end
        end

        function obj = mechanical_parameters(obj, rho, l, m)
            obj.Density = rho;
            obj.lambda = l;
            obj.mu     = m;
            for p=1:obj.numTriangles
                obj.Triangles(p).Density = rho;
                obj.Triangles(p).lambda  = l;
                obj.Triangles(p).mu      = m;
            end
            for p=1:obj.numRectangles
                obj.Rectangles(p).Density = rho;
                obj.Rectangles(p).lambda  = l;
                obj.Rectangles(p).mu      = m;
            end
        end
        function obj = viscous_parameters(obj, lv, mv)
            obj.lambda_vis = lv;
            obj.mu_vis     = mv;
            for p=1:obj.numTriangles
                obj.Triangles(p).lambda_vis  = lv;
                obj.Triangles(p).mu_vis      = mv;
            end
            for p=1:obj.numRectangles
                obj.Rectangles(p).lambda_vis  = lv;
                obj.Rectangles(p).mu_vis      = mv;
            end
        end
        function energy = total_strain_potential_energy(obj, disps)
            energy = 0;
            for p=1:obj.numTriangles
               tri = obj.Triangles(p);
               energy = energy + tri.partial_strain_potential_energy(disps);
            end
            for p=1:obj.numRectangles
               rct = obj.Rectangles(p);
               energy = energy + rct.partial_strain_potential_energy(disps);
            end
        end

        function energy = total_strain_potential_energy_Green_strain(obj, disps)
            energy = 0;
            for p=1:obj.numTriangles
               tri = obj.Triangles(p);
               energy = energy + tri.partial_strain_potential_energy_Green_strain(disps);
            end
            for p=1:obj.numRectangles
               rct = obj.Rectangles(p);
               energy = energy + rct.partial_strain_potential_energy_Green_strain(disps);
            end
        end

        function energy = total_gravitational_potential_energy(obj, disps, grav)
            energy = 0;
            for p=1:obj.numTriangles
               tri = obj.Triangles(p);
               energy = energy + tri.partial_gravitational_potential_energy(disps, grav);
            end
        end

        function output = draw( obj, disps, color )
            arguments
                obj;
                disps = zeros(2, obj.numNodalPoints);
                color = [0.9, 0.9, 0.9];
            end
            for p=1:obj.numTriangles
                tri = obj.Triangles(p);
                vs = tri.Vertices;
                i = vs(1); j = vs(2); k = vs(3);
                np = obj.NodalPoints;
                xi = np(i).Coordinates + disps(:,i);
                xj = np(j).Coordinates + disps(:,j);
                xk = np(k).Coordinates + disps(:,k);
                points = [ xi, xj, xk, xi ];
                plot(points(1,:), points(2,:), 'k-');
                fill(points(1,:), points(2,:), color);
                hold on;
            end
        end

        function output = draw_individual( obj, disps, color )
            arguments
                obj;
                disps = zeros(2, obj.numNodalPoints);
                color = [0.9, 0.9, 0.9];
            end
            for p=1:obj.numTriangles
                tri = obj.Triangles(p);
                vs = tri.Vertices;
                i = vs(1); j = vs(2); k = vs(3);
                np = obj.NodalPoints;
                xi = np(i).Coordinates + disps(:,i);
                xj = np(j).Coordinates + disps(:,j);
                xk = np(k).Coordinates + disps(:,k);
                points = [ xi, xj, xk, xi ];
                plot(points(1,:), points(2,:), 'k-');
                if isempty(obj.Triangles(p).color)
                    fill(points(1,:), points(2,:), color);
                else
                    fill(points(1,:), points(2,:), obj.Triangles(p).color)
                end
                hold on;
            end
        end

        function output = draw_individual_rectangles( obj, disps, color )
            arguments
                obj;
                disps = zeros(2, obj.numNodalPoints);
                color = [0.9, 0.9, 0.9];
            end
            for p=1:obj.numRectangles
                rct = obj.Rectangles(p);
                vs = rct.Vertices;
                i = vs(1); j = vs(2); k = vs(3); l = vs(4);
                np = obj.NodalPoints;
                xi = np(i).Coordinates + disps(:,i);
                xj = np(j).Coordinates + disps(:,j);
                xk = np(k).Coordinates + disps(:,k);
                xl = np(l).Coordinates + disps(:,l);
                points = [ xi, xj, xk, xl, xi ];
                plot(points(1,:), points(2,:), 'k-');
                if isempty(obj.Rectangles(p).color)
                    fill(points(1,:), points(2,:), color);
                else
                    fill(points(1,:), points(2,:), obj.Rectangles(p).color)
                end
                hold on;
            end
        end
        
        function obj = calculate_stiffness_matrix(obj)
            obj.Stiffness_Matrix = zeros(2*obj.numNodalPoints, 2*obj.numNodalPoints);
            %
            for p=1:obj.numTriangles
                tri = obj.Triangles(p);
                vs = tri.Vertices;
                i = vs(1); j = vs(2); k = vs(3);
                [obj.Triangles(p), K_p] = tri.partial_stiffness_matrix;
                loc = [ 2*i-1, 2*i, 2*j-1, 2*j, 2*k-1, 2*k ];
                obj.Stiffness_Matrix(loc,loc) = obj.Stiffness_Matrix(loc,loc) + K_p;
            end
            %
            for p=1:obj.numRectangles
                rct = obj.Rectangles(p);
                vs = rct.Vertices;
                i = vs(1); j = vs(2); k = vs(3); l = vs(4);
                [obj.Rectangles(p), K_p] = rct.partial_stiffness_matrix;
                loc = [ 2*i-1, 2*i, 2*j-1, 2*j, 2*k-1, 2*k, 2*l-1, 2*l ];
                obj.Stiffness_Matrix(loc,loc) = obj.Stiffness_Matrix(loc,loc) + K_p;
            end
        end

        function obj = calculate_damping_matrix(obj)
            obj.Damping_Matrix = zeros(2*obj.numNodalPoints, 2*obj.numNodalPoints);
            %
            for p=1:obj.numTriangles
                tri = obj.Triangles(p);
                vs = tri.Vertices;
                i = vs(1); j = vs(2); k = vs(3);
                [obj.Triangles(p), B_p] = tri.partial_damping_matrix;
                loc = [ 2*i-1, 2*i, 2*j-1, 2*j, 2*k-1, 2*k ];
                obj.Damping_Matrix(loc,loc) = obj.Damping_Matrix(loc,loc) + B_p;
            end
            %
            for p=1:obj.numRectangles
                rct = obj.Rectangles(p);
                vs = rct.Vertices;
                i = vs(1); j = vs(2); k = vs(3); l = vs(4);
                [obj.Rectangles(p), B_p] = rct.partial_damping_matrix;
                loc = [ 2*i-1, 2*i, 2*j-1, 2*j, 2*k-1, 2*k, 2*l-1, 2*l ];
                obj.Damping_Matrix(loc,loc) = obj.Damping_Matrix(loc,loc) + B_p;
            end
        end
        
        function obj = calculate_inertia_matrix(obj)
            obj.Inertia_Matrix = zeros(2*obj.numNodalPoints, 2*obj.numNodalPoints);
            %
            for p=1:obj.numTriangles
                tri = obj.Triangles(p);
                vs = tri.Vertices;
                i = vs(1); j = vs(2); k = vs(3);
                [obj.Triangles(p), M_p] = tri.partial_inertia_matrix;
                loc = [ 2*i-1, 2*i, 2*j-1, 2*j, 2*k-1, 2*k ];
                obj.Inertia_Matrix(loc,loc) = obj.Inertia_Matrix(loc,loc) + M_p;
            end
            %
            for p=1:obj.numRectangles
                rct = obj.Rectangles(p);
                vs = rct.Vertices;
                i = vs(1); j = vs(2); k = vs(3); l = vs(4);
                [obj.Rectangles(p), M_p] = rct.partial_inertia_matrix;
                loc = [ 2*i-1, 2*i, 2*j-1, 2*j, 2*k-1, 2*k, 2*l-1, 2*l ];
                obj.Inertia_Matrix(loc,loc) = obj.Inertia_Matrix(loc,loc) + M_p;
            end
        end

        function obj = calculate_gravitational_vector(obj, g)
            obj.Gravitational_Vector = zeros(2*obj.numNodalPoints, 1);
            for p=1:obj.numTriangles
                tri = obj.Triangles(p);
                vs = tri.Vertices;
                i = vs(1); j = vs(2); k = vs(3);
                [obj.Triangles(p), grav_p] = tri.partial_gravitational_vector(g);
                loc = [ 2*i-1, 2*i, 2*j-1, 2*j, 2*k-1, 2*k ];
                obj.Gravitational_Vector(loc) = obj.Gravitational_Vector(loc) + grav_p;
            end
        end

        function A = constraint_matrix(obj, index)
            index = index(index <= obj.numNodalPoints);
            A = zeros(2*obj.numNodalPoints, 2*length(index));
            for k=1:length(index)
                A(2*index(k)-1, 2*k-1) = 1;
                A(2*index(k),   2*k)   = 1;
            end
        end
        
        function [index, index_rects, index_npoints] = extract_index_for_subregion(obj, index, index_rects)
            for r=1:obj.numSubRegions
                common = intersect(index, obj.SubRegions(r).Index_Triangles);
                if common
                    fprintf("%d ", common);
                    fprintf('already in another subregion\n');
                    index = setdiff(index, common);
                end
                %
                common_rects = intersect(index_rects, obj.SubRegions(r).Index_Rectangles);
                if common_rects
                    fprintf("%d ", common_rects);
                    fprintf('rectangles already in another subregion\n');
                    index_rects = setdiff(index_rects, common_rects);
                end
            end
            index_npoints = [];
            if index
                for p=index
                    tri = obj.Triangles(p);
                    vs = tri.Vertices;
                    index_npoints = union(index_npoints, vs);
                end
            end
            if index_rects
                for p=index_rects
                    rect = obj.Rectangles(p);
                    vs = rect.Vertices;
                    index_npoints = union(index_npoints, vs);
                end
            end
        end

        function obj = define_subregion(obj, index, index_rects)
            arguments
                obj; index;
                index_rects = [];
            end

            if isempty(obj.numSubRegions)
                obj.numSubRegions = 0;
            end
            [index, index_rects, index_npoints] = obj.extract_index_for_subregion(index, index_rects);
            
            obj.numSubRegions = obj.numSubRegions + 1;
            obj.SubRegions = [ obj.SubRegions, SubRegion(index, index_rects, index_npoints) ];
            obj = obj.subregion_partial_connection_matrices;
        end

        function obj = subregion_partial_connection_matrices(obj)
            r = obj.numSubRegions;
            npoints = obj.SubRegions(r).numNodalPoints;
            pJl = zeros(2*npoints, 2*npoints);
            pJm = zeros(2*npoints, 2*npoints);
            
            list = obj.SubRegions(r).Index_NodalPoints;
            %%obj.SubRegions(r).Index_Triangles
            for p=obj.SubRegions(r).Index_Triangles
                tri = obj.Triangles(p);
                vs = tri.Vertices;
                i = vs(1); j = vs(2); k = vs(3);
                ip = find(list==i); jp = find(list==j); kp = find(list==k);
                loc = [ 2*ip-1, 2*ip, 2*jp-1, 2*jp, 2*kp-1, 2*kp ];
                pJl(loc,loc) = pJl(loc,loc) + tri.Partial_J_lambda;
                pJm(loc,loc) = pJm(loc,loc) + tri.Partial_J_mu;
            end

            obj.SubRegions(r).Partial_J_lambda = pJl;
            obj.SubRegions(r).Partial_J_mu     = pJm;
        end

        function obj = subregion_mechanical_parameters(obj, rho, l, m, k)
            arguments
                obj; rho; l; m;
                k = obj.numSubRegions;
            end

            obj.SubRegions(k).Density = rho;
            obj.SubRegions(k).lambda  = l;
            obj.SubRegions(k).mu      = m;
            for p = obj.SubRegions(k).Index_Triangles
                obj.Triangles(p).Density = rho;
                obj.Triangles(p).lambda  = l;
                obj.Triangles(p).mu      = m;
            end
        end

        function obj = subregion_viscous_parameters(obj, lv, mv, k)
            arguments
                obj; lv; mv;
                k = obj.numSubRegions;
            end
            obj.SubRegions(k).lambda_vis  = lv;
            obj.SubRegions(k).mu_vis      = mv;
            for p = obj.SubRegions(k).Index_Triangles
                obj.Triangles(p).lambda_vis  = lv;
                obj.Triangles(p).mu_vis      = mv;
            end
        end

        function obj = subregion_color(obj, cl, k)
            arguments
                obj; cl;
                k = obj.numSubRegions;
            end
            obj.SubRegions(k).color = cl;
            if obj.SubRegions(k).Index_Triangles
                for p = obj.SubRegions(k).Index_Triangles
                    obj.Triangles(p).color = cl;
                end
            end
            if obj.SubRegions(k).Index_Rectangles
                for p = obj.SubRegions(k).Index_Rectangles
                    obj.Rectangles(p).color = cl;
                end
            end
        end

        function forces = nodal_forces_Cauchy_strain_original(obj, disps)
            forces = zeros(2*obj.numNodalPoints,1);
            for p=1:obj.numTriangles
                tri = obj.Triangles(p);
                vs = tri.Vertices;
                i = vs(1); j = vs(2); k = vs(3);
                [fi, fj, fk] = tri.nodal_forces_Cauchy_strain(disps(:,i), disps(:,j), disps(:,k));
                forces(2*i-1:2*i) = forces(2*i-1:2*i) + fi;
                forces(2*j-1:2*j) = forces(2*j-1:2*j) + fj;
                forces(2*k-1:2*k) = forces(2*k-1:2*k) + fk;
            end
        end

        function forces = nodal_forces_Cauchy_strain(obj, disps)
            un = reshape(disps, [2*obj.numNodalPoints, 1]);
            forces = - obj.Stiffness_Matrix*un;
        end

        function forces = nodal_forces_Green_strain_original(obj, disps)
            forces = zeros(2*obj.numNodalPoints,1);
            for p=1:obj.numTriangles
                tri = obj.Triangles(p);
                vs = tri.Vertices;
                i = vs(1); j = vs(2); k = vs(3);
                [fi, fj, fk] = tri.nodal_forces_Green_strain(disps(:,i), disps(:,j), disps(:,k));
                forces(2*i-1:2*i) = forces(2*i-1:2*i) + fi;
                forces(2*j-1:2*j) = forces(2*j-1:2*j) + fj;
                forces(2*k-1:2*k) = forces(2*k-1:2*k) + fk;
            end
        end
        
        function obj = calculate_coefficient_matrices_for_Green_strain(obj)
            [ obj.coef_order_one, obj.coef_order_two, obj.c2i, obj.c2j, obj.coef_order_three, obj.c3i, obj.c3j, obj.c3k ] = ...
                calcuclate_coefficient_matrices_for_Green_strain_forces (obj);
        end
	
        
        
        