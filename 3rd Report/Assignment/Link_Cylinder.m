classdef Link_Cylinder < Link
    properties
        radius; density;
    end
    methods
        function obj = Link_Cylinder (l, r, d)
            obj@Link(l, l/2, 0, 0, 0);
            obj.radius = r;
            obj.density = d;
            obj.mass = d * l * (pi*r^2);
            obj.inertia_of_moment_center = (1/12) * obj.mass * (3*r^2 + l^2);
            obj.inertia_of_moment = obj.inertia_of_moment_center ...
                                  + obj.mass * (obj.length - obj.length_center)^2;
        end
    end
end