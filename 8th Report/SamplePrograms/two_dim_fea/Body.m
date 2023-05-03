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