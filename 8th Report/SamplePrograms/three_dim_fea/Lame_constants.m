function [ lambda, mu ] = Lame_constants ( E, nu )
% calculating Lame's constants
% ƒ‰ƒ‚Ì’è”‚ğŒvZ
    lambda = nu*E/(1+nu)/(1-2*nu);
    mu = E/2/(1+nu);
end
