function bifurcation_information = IK_run_ab_coordinates(I)
% run simulation of Izhikevich model for given current values (or value) to
% find the coordinate of (a, b) where bifurcation happens. The coordinates
% of (a, b) were searched within the given range in Izhikevich (2003)
% Fig.2. where 0 < a < .2 and 0 < b < .3.
% 
% Inputs:
%   - I:    current value (array)
% 
% Outputs:
%   - bifurcation_information: stores the bifurcation information in
%                              structure(object) type. All information are
%                              stored in a same order (index) within in
%                              array
%           |
%           |----------.fxpt_type_current: cell of fixed point types (for current iteration)
%           |----------.fxpt_type_previous: cell of fixed point types (for previous iteration)
%           |----------.a_b: coordinate (a, b)
%           |----------.tau: tau value array
%           |----------.delta: delta value array
%           |----------.I: external current
%           |----------.bifurcation_type: type of bifurcation

RelTol = 1e-10; c = -65; d = 2;
A = 0:.001:.2; % array of parameter 'a' of Izhikevich model
B = 0:.001:.3; % array of parameter 'b' of Izhikevich model

num_fxpt = 0; fxpt_type_old = {'dummy','dummy2'};
bifurcation_information = [];

for i = 1:numel(I)   
    % find (a, b) where bifurcation occur    
    for a = A
        j= 1;
        for b = B
            pars = [a, b, c, d, I(i)];
            fxpt_type_current = {};
            delta_current = []; tau_current = [];           
            
            rts = roots_(pars); % v* value of fixed points
            rts_real = find(abs(imag(rts)) < RelTol); % real fxpt index
            [~, unique_idx] = unique(rts); % check if only one fxpt exist
            rts_real = intersect(rts_real, unique_idx)';
            
            for idx = rts_real
                v = rts(idx); u = b*v; % fxpt coordinates
                A = Izhikevich_Jacobian(v, u, pars); % Jacobian of Izhikevich model
                [fxpt_type, ~, ~, delta, tau] = type_(A);
                
                fxpt_type_current{end+1} = fxpt_type; % store fxpt types
                delta_current(end+1) = delta; tau_current(end+1) = tau;
                
            end
            
            % check if bifurcation happens
            % cases:
            %   1) when fxpt number change
            %   2) if fxpt type change            
            if (num_fxpt ~= numel(rts_real)) && (j ~= 1) % case 1)
                
                bifurcation_information(end+1).fxpt_type_current = fxpt_type_current;
                bifurcation_information(end).fxpt_type_previous = fxpt_type_old;
                bifurcation_information(end).a_b = [a, b];
                bifurcation_information(end).tau = tau_current;
                bifurcation_information(end).delta = delta_current;
                bifurcation_information(end).I = I(i);
                bifurcation_information(end).bifurcation_type = 'change in number of fxpt';
                
            elseif ~isempty(setdiff(fxpt_type_old, fxpt_type_current)) && (j ~= 1) % case 2)
                
                bifurcation_information(end+1).fxpt_type_current = fxpt_type_current;
                bifurcation_information(end).fxpt_type_previous = fxpt_type_old;
                bifurcation_information(end).a_b = [a, b];
                bifurcation_information(end).tau = tau_current;
                bifurcation_information(end).delta = delta_current;
                bifurcation_information(end).I = I(i);
                bifurcation_information(end).bifurcation_type = 'change in type of fxpt';
                
            end
            
            num_fxpt = numel(rts_real); fxpt_type_old = fxpt_type_current;
            j = j + 1;         
        end
    end
    
end

end