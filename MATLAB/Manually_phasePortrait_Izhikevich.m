close all;clc; clear; % clear cmd window and cache
addpath('./func') % add path of funtions
%%  basic settings
a = .2958; b = .263; c = -65; d = 2; RelTol = 1e-10;
I0 = 0.2;   % potentiation (steady state) 
I1 = 0.4;  % turn on
I2 = 0;     % turn off
pars = [a b c d I0];

%   phase plane ranges
v1=-80; v2=40;    dv=(v2-v1)/10;
u1=-18; u2=-10;    du=(u2-u1)/10;
v = v1 : dv/100 : v2;
%%   Get trajectories
I_list = [0.2 0.4 0.2 0 0.2 0];
T_stim = 200;
T_list = [100 200 250+T_stim 450+T_stim 500+T_stim 1000+T_stim];
% I_list = [0.4 0.4];
% T_stim = 200;
% T_list = [0.1 400];
% varargin = {'tspan', 1000, 'delta', .01, 'a', .1, 'b', .26, 'c', -65, ...
%         'd', 2, 'I', I_list, 'injectionTime', T_list};
varargin = {'tspan', 1000, 'delta', .01, 'a', pars(1), 'b', pars(2), 'c', pars(3), ...
        'd', pars(4), 'I', I_list, 'injectionTime', T_list};
[tout, xout, teout, event_type] = Izhikevich(-65, -16, varargin{:});
figure(1);
yyaxis left;
plot(tout(2:end-1),xout(2:end-1,1));
ylabel('v');
yyaxis right;
stairs([0,T_list],[0,I_list]);
ylabel('I');
set(gcf,'position',[100,100,1000,250]);
xlabel('Time (ms)');
%%  Phase portrait for different stages
I_list_ext = [0,I_list];
T_list_ext = [0,T_list];
for i = 1: length(I_list)
    figure(i+1);box on; hold on;
    pars(5) = I_list_ext(i);
    [vnull, unull] = Izhikevich_null(v, pars);
    plot(v,vnull,'r-',v,unull,'-g');
    xlabel('u'); ylabel('v');
%     title(['u-v phase plane of IK model ' num2str(pars)]);
    % 	vector field
    [V,U] = meshgrid(v1:dv:v2, u1:du:u2);
    Zv = 0.04*V.^2+5*V+140-U+pars(5);
    Zu = pars(1)*(pars(2)*V-U);
    
    %   normalized quiver
    arrow_len = 1;
    Z_norm_v = Zv./sqrt(Zv.^2+Zu.^2);
    Z_norm_u = Zu./sqrt(Zv.^2+Zu.^2);
    quiver(V,U,Z_norm_v*arrow_len,Z_norm_u*arrow_len,'AutoScale','off');
    %   scaled quiver
%     quiver(V,U,Zv*0.015,Zu*0.015,'AutoScale','off');
    %   roots & types
    Vrts = roots_(pars);
    Urts = (pars(2)*Vrts);
    disp(['---------' num2str(i) '--------']);
    for k = 1:length(Vrts)
        J = Izhikevich_Jacobian(Vrts(k), Urts(k), pars);
        [str, ~, d, delta, tau] = type_(J);
        
        if strcmp(str(1:6),'stable')
            scatter(Vrts(k),Urts(k),'ko','filled');
        else
            scatter(Vrts(k),Urts(k),'ro','filled');
        end
        disp([str,' ', num2str(Vrts(k))]);
    end
    idx_sel = find(T_list_ext(i)<tout & tout<T_list_ext(i+1));
    plot(xout(idx_sel,1),xout(idx_sel,2),'k-');
    xlim([v1, v2]);ylim([u1, u2]);
    set(gcf,'position',[100,100,400,300]);
end

