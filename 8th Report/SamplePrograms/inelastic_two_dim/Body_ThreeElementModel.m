classdef Body_ThreeElementModel < Body
    properties
        lambda_vis_1, mu_vis_1, lambda_vis_2, mu_vis_2;
    end

    methods
        function obj = Body_ThreeElementModel(npoints, points, ntris, tris, h)
            obj@Body(npoints, points, ntris, tris, h);
            for p=1:ntris
                i = tris(p,1); j = tris(p,2); k = tris(p,3);
                tr(p) = Triangle_ThreeElementModel(i, j, k, points(:,i), points(:,j), points(:,k), h);
            end
            obj.Triangles = tr;
        end