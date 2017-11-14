function [dvdt, dudt]  = Izhikevich_VF(v, u, pars)
% vector field of HR model
% Inputs:
%   - v: membrane potential variable
%   - u: recovery varaible
%   - pars: parameters for Izhikevich model
%           [a b c d I]

dvdt = 0.04.*v.^2 + 5.*v + 140 - u + pars(5);
dudt = pars(1).*(pars(2).*v - u);

end