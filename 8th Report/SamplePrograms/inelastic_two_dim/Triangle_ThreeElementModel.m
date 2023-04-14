classdef Triangle_ThreeElementModel < Triangle
    properties
        lambda_vis_1, mu_vis_1, lambda_vis_2, mu_vis_2;
    end
    methods
        function obj = mechanical_parameters(obj, rho, l, m, lv1, mv1, lv2, mv2)
            obj = obj.mechanical_parameters@Triangle(rho, l, m);
            obj.lambda_vis_1 = lv1; obj.mu_vis_1 = mv1;
            obj.lambda_vis_2 = lv2; obj.mu_vis_2 = mv2;
        end
    end
end
