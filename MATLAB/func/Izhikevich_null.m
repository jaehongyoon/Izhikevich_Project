function [vnull, unull] = Izhikevich_null(v, pars)
% nullclines of Izhikevich model

vnull = .04.*v.^2+5.*v+140-pars(5); % u = vnull
unull = pars(2).*v; % u = unull

end