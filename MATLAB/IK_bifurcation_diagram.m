clc; clear; % clear cmd window and cache
addpath('./func') % add path of funtions

%% ==================================
% Analyze the fixed point for given parameter sets
% ===================================

I = .3; c = -65; d = 2; RelTol = 1e-10; % THIS PART SHOULD CHANGE ACCORDINGLY 
                                        % TO THE SPECIFICATION OF ANALYSIS
                                        
% structure (object) of information storage
% fxpt
% |
% |------.a: value of a
% |------.b: value of b
% |------.type: either 1 or 2
% |------.v_star: v* value array
% |------.I: external current I
% |------.bf_type: bifucation type

fxpt = [];
da = .001; A = 0:da:.15;
db = .001; B = .1:db:.3;

h = waitbar(0, 'bifurcation calculating');
for a = A
    waitbar(a/max(A));
    for b = B
        pars = [a, b, c, d, I]; % parameters for the simulation
        rts = roots_(pars); % v* value of fixed points
        rts_real = find(abs(imag(rts)) < RelTol); % regard those whose imaginary part is less than tolerance as real fixed point
        
        [~, unique_idx] = unique(rts); % select only unique fxpt points
        rts_real = intersect(rts_real, unique_idx)';
        
        for idx = rts_real
            v = rts(idx); u = b*v; % fxpt coordinate
            J = Izhikevich_Jacobian(v, u, pars); % Jacobian of Izhikevich model
            [bf_type, ~, ~, ~, ~] = type_(J);  % fxpt types
        
            % check if the current data type is in fxpt array
            [TF, j] = in_structure2(fxpt, idx, bf_type, a, I);

            if TF % if inside array, add x*/I values to fxpt
                fxpt(j).v_star(end+1) = v; fxpt(j).b(end+1) = b;         
            else % if not exist, add to fxpt array
                fxpt(end+1).type = idx; fxpt(end).bf_type = bf_type;
                fxpt(end).v_star = v; fxpt(end).I = I;
                fxpt(end).a = a; fxpt(end).b = b;
            end
        end        
    end
end

close(h);

%% ==================================
% Analyze the fixed point for given parameter sets
% ===================================

filename = 'result.gif';
sorted = bifurcation_animation(fxpt, filename);



