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