clc; clear; % clear cmd window and cache
addpath('../func') % add path of funtions

%% ==================================
% Analyze the fixed point for given parameter sets
% ===================================

a = .1; b = .2; c = -65; d = 2; RelTol = 1e-10; % THIS PART SHOULD CHANGE 
                                                % ACCORDINGLY TO THE
                                                % SPECIFICATION OF ANALYSIS

% structure (object) of information storage
% fxpt
% |
% |------.type: either 1 or 2
% |------.v_star: v* value array
% |------.pars: parameters of Izhikevich model including external current I
% |------.bf_type: bifucation type
% |------.delta: delta of Jacobian A
% |------.tau: tau of Jacobian A

% use structure (object) for storing information of trajectory
% trajectory 
% |
% |----------.fxpt_type: either 1 or 2
% |----------.pars: parameters of Izhikevich model including external current I
% |----------.evector: eigenvector of given fixed point
% |----------.evalue: eigenvalue of given fixed point
% |----------.v1: membrane potential variable trajectory near first fxpt
% |----------.u1: recovery variable trajectory near first txpt
% |----------.v2: membrane potential variable trajectory near second fxpt
% |----------.u2: recovery variable trajectory near second txpt

fxpt = []; trajectory = [];
h = waitbar(0, 'bifurcation calculating');
dI = .01; I = -1.5:dI:3;
tspan = [0 1000];

for i = 1:numel(I)
    waitbar(i/numel(I));
    
    pars = [a, b, c, d, I(i)]; % parameters for current simulation
    paramset = {'tspan', tspan, 'delta', .01, 'a', a, 'b', b, 'c', c, ...
        'd', d, 'I', I(i), 'injectionTime', [0]};
    
    trajectory(i).pars = pars; % store parameters for current simulation to trajectory structure
    % initialize trajectory object for given iteration
    trajectory(i).evector = {'NaN', 'NaN'};
    trajectory(i).evalue = {'NaN', 'NaN'};
    
    rts = roots_(pars); % v* value of fixed points
    rts_real = find(abs(imag(rts)) < RelTol); % real fxpt index
    
    % check if only one fxpt exist
    [~, unique_idx] = unique(rts);
    rts_real = intersect(rts_real, unique_idx)';
    
    % fixed point analysis
    for idx = rts_real
        if isfield(trajectory(i), 'fxpt_type')
            trajectory(i).fxpt_type(end+1) = idx;
        else
            trajectory(i).fxpt_type = idx;
        end
        
        v = rts(idx); u = b*v; % fxpt coordinate
        
        % fxpt types
        A = Izhikevich_Jacobian(v, u); % Jacobian of Izhikevich model
        [bf_type, eigenvector, eigenvalue, delta, tau] = type_(A);
        
        % retrieve trajectory information
        trajectory(i).evector{idx} = eigenvector; % eigenvector
        trajectory(i).evalue{idx} = eigenvalue; % eigenvalue
        
        % check if the current data type is in fxpt array
        [TF, i] = in_structure(fxpt, idx, bf_type);
        
        if TF % if inside array, add x*/I values to fxpt
            fxpt(i).v_star(end+1) = v; fxpt(i).pars{end+1} = pars;
            fxpt(i).delta(end+1) = delta; fxpt(i).tau(end+1) = tau;            
        else % if not exist, add to fxpt array
            fxpt(end+1).type = idx; fxpt(end).bf_type = bf_type;
            fxpt(end).v_star = v; fxpt(end).pars = {pars};
            fxpt(end).delta = delta; fxpt(end).tau = tau;
        end
        
        % get trajectory near fixed point
        [tout, xout, ~, ~] = Izhikevich(v, u+.1, paramset);
        if idx == 1
            trajectory(i).v1 = xout(:,1);
            trajectory(i).u1 = xout(:,2);
        elseif idx == 2
            trajectory(i).v2 = xout(:,1);
            trajectory(i).u2 = xout(:,2);
        end
    end
    
end

close(h);

%% ==================================
% plot the results
%   bifurcation graph, tau-delta graph, phase diagram
% ===================================

Izhikevich_Results_Plot(fxpt, trajectory);

return