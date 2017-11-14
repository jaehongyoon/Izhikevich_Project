function [tout, xout, teout, event_type] = Izhikevich(vo, uo, varargin)
% simulation of Izhikevich model
% Inputs:
%   - vo: initial value of membrane voltage variable
%   - u0: initial value of recovery variable
% 
% varargin:
%   - 'tspan': time span of simulation in vector form [t start, t end]
%   - 'delta': dt of simulation
%   - 'a', 'b', 'c', 'd': 4 varaibles of Izhikevich model 
%   - 'I': external current in a vector format [I0, I1, I2, ...]
%   - 'injectionTime': time when the currents were injected
%                      if current injection end before tspan, the user
%                      have to add the injection_end time as well with
%                      I(injection_end) = 0
% Outputs:
%   - tout: simulation time array
%   - xout: simulation result of v(t) and u(t)
%   - teout: event time log
%   - event_type: type of event happend in recorded time

% parse function inputs
funcInputs = parseMyInputs(varargin{:});

% parameter array
pars = [funcInputs.a, funcInputs.b, funcInputs.c, funcInputs.d, 0, ];

% initial conditions
xo_ = [vo, uo];

if isempty(funcInputs.injectionTime)
    tend = funcInputs.tspan(2); % if no current injection was given, run until the end of the tspan
else
    tend = funcInputs.injectionTime(1); % run without current until first injection was given
end
tspan = funcInputs.tspan(1):funcInputs.delta:tend;
options = odeset('RelTol', 1e-10);
options = odeset(options, 'Event', @(t, y) events(t, y));

% output vectors
tout = funcInputs.tspan(1);
xout = xo_; teout = []; event_type = {}; next = 1;

while(1) % run simulation
    [t_, x_, te_, ~, ~] = ode45(@(t, x) Izhikevich_ODE(t, x, pars), tspan, xo_, options);
    xo_ = x_(end,:); % the last value of simulation
    tout = [tout; t_(2:end-1)]; % tarray of simulation
    xout = [xout; x_(2:end-1, :)]; % v(t), u(t) array
    teout = [teout; te_]; % array of event time
    
    if(te_) % if spike occured
        xo_ = [funcInputs.c, xo_(2)+funcInputs.d, ];
        event_type{end+1} = 'spike';
    else
        if t_(end) < funcInputs.tspan(2) % if the simulation ended before end of timespan
            % new injection
            I = funcInputs.I(next);
            pars = [funcInputs.a, funcInputs.b, funcInputs.c, funcInputs.d, I, ];
            next = next+1;
            event_type{end+1} = 'current injection';
            if next < numel(funcInputs.I) % if current injection changes further more
                tend = funcInputs.injectionTime(next);
            else
                tend = funcInputs.tspan(2);
            end
        else % end of simulation
            break
        end
    end
    
    % resume simulation
    options = odeset(options, 'InitialStep', t_(end)-t_(end-1), 'MaxStep', t_(end)-t_(1));
    tspan = [t_(end), tout(end)+dt:dt:tend];
end
    
end

function rhs = Izhikevich_ODE(t, state, pars)
rhs(1) = 0.04.*state(1).^2 + 5.*state(1) + 140 - state(2) + pars(5);
rhs(2) = pars(1).*(pars(2).*state(1) - state(2));
end

function [value,isterminal,direction] = events(t, state) % ODE event location 
    value = [1,];
    isterminal = [1,];
    direction = [1,];
   
    v = state(1);
    u = state(2);

    % event of ODE = spike
    value(1) = (v > 30) - 1; % when spike

end

function funcInputs = parseMyInputs(varargin)
funcInputs = inputParser;

addRequired(funcInputs, 'tspan', [0 1000], @isvector);
addRequired(funcInputs, 'delta', .01, @isscalar);
addParameter(funcInputs, 'a', .1, @isscalar);
addParameter(funcInputs, 'b', .2, @isscalar);
addParameter(funcInputs, 'c', -65.0, @isscalar);
addParameter(funcInputs, 'd', 2.0, @isscalar);
addParameter(funcInputs, 'I', [], @isvector);
addParameter(funcInputs, 'injectionTime', [], @isvector);

parse(funcInputs, varargin{:});
funcInputs = funcInputs.Results;
end