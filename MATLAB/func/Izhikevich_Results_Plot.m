function Izhikevich_Results_Plot(fxpt, trajectory)

%% ==========================================
% Bifurcation Diagram
% ===========================================
fig1 = figure(100);
hold on
for i = 1:length(fxpt)
    bf_type = fxpt(i).bf_type;
    if (~isempty(strfind(bf_type, 'stable'))) &&...
            (isempty(strfind(bf_type, 'unstable'))) % if a stable fxpt
        line_type = '-';
    else % if unstable, use dotted line
        line_type = '--';
    end
    
    if fxpt(i).type == 1
        color_spec = 'b';
    elseif fxpt(i).type == 2
        color_spec = 'r';
    end
    
    str = sprintf('%s, %s', fxpt(i).type, fxpt(i).bf_type);
    I = zeros(1, length(fxpt(i).pars));
    for j = 1:length(fxpt(i).pars)
        I(j) = fxpt(i).pars{j}(5);
    end
    
    if length(I) > 1
        plot(I, fxpt(i).v_star, [color_spec, line_type], ...
            'linewidth', 0.5) % draw line
    else
        plot(fxpt(i).I, fxpt(i).x_star, [color_spec, '*'], ...
            'linewidth', 0.5) % draw line
    end
end
xlabel('current I'); ylabel('v* of fixed point')
hold off

%% ==========================================
% Delta-tau graph
% ===========================================
fig2 = figure(200);
hold on 
delta = 0:0.01:4.5;
plot(delta, 2*sqrt(delta), 'k-')
plot(delta, -2*sqrt(delta), 'k-')
for i = 1:length(fxpt)    
    bf_type = fxpt(i).bf_type;    
    if (~isempty(strfind(bf_type, 'stable'))) &&...
            (isempty(strfind(bf_type, 'unstable'))) % if a stable fxpt
        line_type = '-';
    else
        line_type = '--';
    end
    
    if fxpt(i).type == 1
        color_spec = 'b';
    elseif fxpt(i).type == 2
        color_spec = 'r';
    end    
           
    eval(sprintf('line_%s', strcat(num2str(fxpt(i).type)),...
        '= plot(fxpt(i).delta, fxpt(i).tau, [color_spec, line_type]);'));    
end
xlabel('delta');ylabel('tau')
hold off

%% ==========================================
% Phase diagram
% ===========================================

dv = 0.01; du = 0.01;
v1 = -100; v2 = 100; u1 = -100; u2 = 100;
v = v1:dv:v2; u = u1:du:u2;
[V,U] = meshgrid(v1:1:v2, u1:1:u2);

fig3 = figure(300);
xlabel('x');ylabel('y')
hold on
for i = 1:length(trajectory)    
    I = trajectory(i).pars(5); % get current value    
    
    [vnull, unull] = Izhikevich_null(v, trajectory(i).pars); % nullclines
    [dvdt, dudt]  = Izhikevich_VF(V, U, trajectory(i).pars); % vector fields
    
    % plot trajectories
    trajectory_1 = plot(trajectory(i).v1, trajectory(i).u1, 'b-');
    trajectory_2 = plot(trajectory(i).v2, trajectory(i).u2, 'r-');
    
    vectorfield = quiver(V, U, dvdt, dudt, .5, 'g'); % plot vector field    
    vnullcline = plot(v, vnull, 'k-');
    unullcline = plot(v, unull, 'k--'); % plot nullclines
        
    pause(); % pause to make decision/observe the change in trajectory
    axis([-100 50 -100 100])
    % delete plots
    delete(vnullcline); delete(unullcline); delete(vectorfield);
    delete(trajectory_1); delete(trajectory_2);

end
hold off

end