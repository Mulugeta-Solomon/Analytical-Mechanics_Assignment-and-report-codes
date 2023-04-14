classdef Body
    properties
        numNodalPoints; NodalPoints;
        numTetrahedra; Tetrahedra;
        BoundaryFaces;
        numSubRegions; SubRegions;
        numSurfaces; Surfaces;
        Density; lambda; mu; lambda_vis; mu_vis;
        strain_potential_energy;
        gravitational_potential_energy;
        J_lambda; J_mu;
        Stiffness_Matrix;
        Damping_Matrix;
        Inertia_Matrix;
        Gravitational_Vector;
    end

    methods
        function obj = Body(npoints, points, ntetras, tetras)
            obj.numNodalPoints = npoints;
            for k=1:npoints
                pt(k) = NodalPoint(points(:,k));
            end
            obj.NodalPoints = pt;
            %
            if ntetras > 0
                obj.numTetrahedra = ntetras;
                for p=1:ntetras
                    i = tetras(p,1); j = tetras(p,2); k = tetras(p,3); l = tetras(p,4);
                    te(p) = Tetrahedron(i, j, k, l, points(:,i), points(:,j), points(:,k), points(:,l));
                end
                obj.Tetrahedra = te;
            end
            obj.J_lambda = zeros(3*obj.numNodalPoints, 3*obj.numNodalPoints);
            obj.J_mu     = zeros(3*obj.numNodalPoints, 3*obj.numNodalPoints);
            for p=1:obj.numTetrahedra
                tetra = obj.Tetrahedra(p);
                vs = tetra.Vertices;
                i = vs(1); j = vs(2); k = vs(3); l = vs(4);
                loc = [ 3*i-2, 3*i-1, 3*i, 3*j-2, 3*j-1, 3*j, 3*k-2, 3*k-1, 3*k, 3*l-2, 3*l-1, 3*l ];
                obj.J_lambda(loc,loc) = obj.J_lambda(loc,loc) + tetra.Partial_J_lambda;
                obj.J_mu(loc,loc)     = obj.J_mu(loc,loc)     + tetra.Partial_J_mu;
            end

            faces = [];
            for p=1:obj.numTetrahedra
                tetra = obj.Tetrahedra(p);
                faces = [ faces; tetra.Faces ];
            end
            faces_sorted = sort(faces, 2);
            n = size(faces_sorted, 1);
            onlyone = 1:n;
            for p=1:n-1
                [lia, loc] = ismember(faces_sorted(p,:), faces_sorted(p+1:n,:), 'rows');
                if lia
                    onlyone(p) = 0;
                    onlyone(p+loc) = 0;
                end
            end
            onlyone = setdiff(onlyone, 0);
            obj.BoundaryFaces = faces(onlyone, :);
            
            [ nsufs, surfs, order ] = triangles_to_surfaces (obj.BoundaryFaces);
            obj.numSurfaces = nsufs;
            for p=1:nsufs
                [surf_p, area_p] = extract_surface(p, surfs, order, obj.BoundaryFaces, points );
                sf(p) = Surface(size(surf_p,1), surf_p, area_p);
            end
            obj.Surfaces = sf;
        end

        function obj = mechanical_parameters(obj, rho, l, m)
            obj.Density = rho;
            obj.lambda = l;
            obj.mu     = m;
            for p=1:obj.numTetrahedra
                obj.Tetrahedra(p).Density = rho;
                obj.Tetrahedra(p).lambda  = l;
                obj.Tetrahedra(p).mu      = m;
            end
        end

        function obj = viscous_parameters(obj, lv, mv)
            obj.lambda_vis = lv;
            obj.mu_vis     = mv;
            for p=1:obj.numTetrahedra
                obj.Tetrahedra(p).lambda_vis  = lv;
                obj.Tetrahedra(p).mu_vis      = mv;
            end
        end

        function energy = total_strain_potential_energy(obj, disps)
            energy = 0;
            for p=1:obj.numTetrahedra
               tetra = obj.Tetrahedra(p);
               energy = energy + tetra.partial_strain_potential_energy(disps);
            end
        end

        function energy = total_strain_potential_energy_Green_strain(obj, disps)
            energy = 0;
            for p=1:obj.numTetrahedra
               tetra = obj.Tetrahedra(p);
               energy = energy + tetra.partial_strain_potential_energy_Green_strain(disps);
            end
        end