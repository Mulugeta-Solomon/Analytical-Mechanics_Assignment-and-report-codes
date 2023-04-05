classdef Link
    properties
        length;
        length_center;
        mass;
        inertia_of_moment_center;
        inertia_of_moment;
    end
    methods
        function obj = Link (l, lc, m, Jc, J)
            obj.length = l;
            obj.length_center = lc;
            obj.mass = m;
            obj.inertia_of_moment_center = Jc;
            obj.inertia_of_moment = J;
        end
    end
end