close all;clc; clear; % clear cmd window and cache
addpath('./func') % add path of funtions
%%  basic settings
a = .1; b = .26; c = -65; d = 2; RelTol = 1e-10;
I0 = 0.2;   % potentiation (steady state)
I1 = 0.35;  % turn on
I2 = 0;     % turn off
%%  limit cycle plot
%   IDEA: try a grid of initial values, see which ones enters limit cycle.
%   then plot boundary between two sets of points
v1=-80; v2=40;    dv=(v2-v1)/10;
u1=-18; u2=-10;   du=(u2-u1)/10;
v = v1 : dv : v2;
u = u1 : du : u2;

[VV,UU] = meshgrid(v,u);
firing = zeros(size(UU));
parfor i = 1:length(firing(:))
    pars = [a b c d I0];
    I_list = [1,1,1]*pars(5);
    T_list = [0.1, 200,300];
    varargin = {'tspan', 300, 'delta', .01, 'a', pars(1), 'b', pars(2), 'c', pars(3), ...
        'd', pars(4), 'I', I_list, 'injectionTime', T_list};
    [tout, xout, teout, event_type] = Izhikevich(VV(i),UU(i), varargin{:});
    % judge if there is sustained firing
    firing(i) = any(xout(tout>200)>=29.7) & any(xout(tout>200)<=-65);
end

figure; hold on;
scatter(VV(firing>0.5),UU(firing>0.5),'r.')
scatter(VV(firing>0.5),UU(firing>0.5),'k.')
xlabel('initial v'); ylabel('initial u');
xlim([v1,v2]);ylim([u1, u2]);
set(gcf,'position',[100,100,400,300]);
