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