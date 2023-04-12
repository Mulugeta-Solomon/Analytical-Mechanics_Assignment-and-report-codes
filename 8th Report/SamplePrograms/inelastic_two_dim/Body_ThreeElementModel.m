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

        function obj = mechanical_parameters(obj, rho, l, m, lv1, mv1, lv2, mv2)
            obj = mechanical_parameters@Body(obj, rho, l, m);
            obj.lambda_vis_1 = lv1; obj.mu_vis_1 = mv1;
            obj.lambda_vis_2 = lv2; obj.mu_vis_2 = mv2;
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
            obj.SubRegions = [ obj.SubRegions, SubRegion_ThreeElementModel(index, index_rects, index_npoints) ];
            obj = obj.subregion_partial_connection_matrices;
            
            obj.SubRegions(obj.numSubRegions).suffixes_for_nodal_points = ...
                reshape( [ 2*index_npoints-1, 2*index_npoints ]', [ 2*length(index_npoints), 1 ] );
        end