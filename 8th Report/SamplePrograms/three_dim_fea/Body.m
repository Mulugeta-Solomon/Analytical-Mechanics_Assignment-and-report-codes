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