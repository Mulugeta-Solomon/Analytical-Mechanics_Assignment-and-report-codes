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

        function obj = subregion_mechanical_parameters(obj, rho, l, m, lv1, mv1, lv2, mv2, k)
            arguments
                obj; rho; l; m, lv1, mv1, lv2, mv2;
                k = obj.numSubRegions;
            end
            obj = obj.subregion_mechanical_parameters@Body(rho, l, m, k);
            obj.SubRegions(k).lambda_vis_1 = lv1; obj.SubRegions(k).mu_vis_1 = mv1;
            obj.SubRegions(k).lambda_vis_2 = lv2; obj.SubRegions(k).mu_vis_2 = mv2;
            for p = obj.SubRegions(k).Index_Triangles
                obj.Triangles(p) = obj.Triangles(p).mechanical_parameters(rho, l, m, lv1, mv1, lv2, mv2);
            end
        end

        function [ps, pe] = flambda_fmu_location_corresponding_to_subregions(obj)
            count = [];
            for p=1:obj.numSubRegions
                count = [ count, obj.SubRegions(p).numNodalPoints ];
            end
            ps(1) = 1; pe(1) = 2*count(1);
            for p=2:obj.numSubRegions
                ps(p) = pe(p-1) + 1;
                pe(p) = ps(p) -1 + 2*count(p);
            end
        end
        