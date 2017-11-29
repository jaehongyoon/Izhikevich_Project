close all;clc; clear; % clear cmd window and cache
addpath('./func') % add path of funtions
%%  basic settings
a = .1; b = .26; c = -65; d = 2; RelTol = 1e-10;
I0 = 0.2;   % potentiation (steady state) 
I1 = 0.35;  % turn on
I2 = 0;     % turn off
pars = [a b c d I0];


v1=-80; v2=40;    dv=(v2-v1)/100;
u1=-18; u2=-10;    du=(u2-u1)/100;
v = v1 : dv/100 : v2;

figure; hold on;

%%  vector field
[V,U] = meshgrid(v1:dv:v2, u1:du:u2);
Zv = 0.04*V.^2+5*V+140-U+pars(5);
Zu = pars(1)*(pars(2)*V-U);
%   normalized quiver
arrow_len = 0.3;
Z_norm_v = Zv./sqrt(Zv.^2+Zu.^2);
Z_norm_u = Zu./sqrt(Zv.^2+Zu.^2);
% quiver(V,U,Z_norm_v*arrow_len,Z_norm_u*arrow_len,'AutoScale','off');

%%  special regions
%   1) downward
idx = find(Zu<0);
k = boundary(V(idx),U(idx));
h = fill(V(idx(k)),U(idx(k)),'g','linestyle','none');set(h,'facealpha',.5);
idx = find(Zv>0);
k = boundary(V(idx),U(idx));
h = fill(V(idx(k)),U(idx(k)),'r','linestyle','none');set(h,'facealpha',.5);

%%  Nullclines 
[vnull, unull] = Izhikevich_null(v, pars);
h_nullcline = plot(v,vnull,'r-',v,unull,'-g');
set(h_nullcline,'linewidth',2);
xlabel('u'); ylabel('v');
%%  Figure settings


xlim([v1, v2]);ylim([u1, u2]);
set(gcf,'position',[100,100,400,300]);
