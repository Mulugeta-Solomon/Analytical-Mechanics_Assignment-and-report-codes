classdef Spring
    properties
        Stiffness;
        Damping;
        Natural_Length;
    end
    methods
        function obj = Spring (k, b, len)
            obj.Stiffness = k;
            obj.Damping = b;
            obj.Natural_Length = len;
        end
        
        function [ fi, fj ] = end_forces (obj, xi, vi, xj, vj, extensional_force)
            dij = xj - xi;
            lij = sqrt(dij'*dij);
            eij = dij/lij;
            v_tangent = ( vj - vi )' * eij;
            fi = ( -extensional_force + obj.Stiffness * (lij - obj.Natural_Length) + obj.Damping * v_tangent )*eij;
            fj = -fi;
        end
        
        function energy = internal_energy(obj, xi, xj, extensional_force)
            dij = xj - xi;
            lij = sqrt(dij'*dij);
            extension = lij - obj.Natural_Length;
            energy = (1/2)*obj.Stiffness*extension^2 - extensional_force*extension;
        end
    end
end
