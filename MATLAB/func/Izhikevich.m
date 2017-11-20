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
    tend = funcInputs.tspan; % if no current injection was given, run until the end of the tspan
else
    tend = funcInputs.injectionTime(1); % run without current until first injection was given
end
tspan = 0:funcInputs.delta:tend;
options = odeset('RelTol', 1e-10);
options = odeset(options, 'Events', @(t, y) events(t, y));

% output vectors
tout = funcInputs.tspan(1);
xout = xo_; teout = []; event_type = {}; next = 1;

while(1) % run simulation
    [t_, x_, te_, xe_, ie_] = ode45(@(t, x) Izhikevich_ODE(t, x, pars), tspan, xo_, options);
    xo_ = x_(end,:); % the last value of simulation
    tout = [tout; t_(2:end-1)]; % tarray of simulation
    xout = [xout; x_(2:end-1, :)]; % v(t), u(t) array
    teout = [teout; te_]; % array of event time
    
    if(te_) % if spike occured
        xo_ = [funcInputs.c, xo_(2)+funcInputs.d, ];
        event_type{end+1} = 'spike';
    else       
        if t_(end) < funcInputs.tspan % if the simulation ended before end of timespan
            % new injection
            event_type{end+1} = 'current injection';
            if next <= numel(funcInputs.I) % if current injection changes further more
                I = funcInputs.I(next);
                if next + 1 <= numel(funcInputs.I)
                    tend = funcInputs.injectionTime(next+1);
                else
                    tend = funcInputs.tspan;
                end
            else
                tend = funcInputs.tspan;
            end
            
            pars = [funcInputs.a, funcInputs.b, funcInputs.c, funcInputs.d, I, ];
            next = next+1;
            
        else % end of simulation
            break
        end
    end
    
    % resume simulation
    options = odeset(options, 'InitialStep', t_(end)-t_(end-1), 'MaxStep', t_(end)-t_(1));
    tspan = [t_(end), tout(end)+funcInputs.delta*1.5:funcInputs.delta:tend];
    
    if length(tspan) == 1
        break
    end
end
    
end

function rhs = Izhikevich_ODE(t, state, pars)
    rhs = zeros(length(state),1);

    v = state(1);
    u = state(2);
    %fprintf('v = %.5f, u = %.5f\n', v, u);

    a = pars(1);
    b = pars(2);
    c = pars(3);
    d = pars(4);
    I = pars(5);

    phi = 0.04 .* v .^ 2 + 5 .* v + 140;

    rhs(1) = phi - u + I;
    rhs(2) = a .* (b .* v - u);
end

function [value,isterminal,direction] = events(t, state) % ODE event location 
    isterminal = 1;
    direction = 1;
   
    v = state(1);
    u = state(2);

    %fprintf('v = %.5f, u = %.5f\n', v, u);
    % event of ODE = spike
    value = (v > 30)-1; % when spike

end

function funcInputs = parseMyInputs(varargin)
funcInputs = inputParser;

addRequired(funcInputs, 'dummy', @checkDummy);
addParameter(funcInputs, 'tspan', 1000, @isscalar);
addParameter(funcInputs, 'delta', .01, @isscalar);
addParameter(funcInputs, 'a', .1, @isscalar);
addParameter(funcInputs, 'b', .2, @isscalar);
addParameter(funcInputs, 'c', -65.0, @isscalar);
addParameter(funcInputs, 'd', 2.0, @isscalar);
addParameter(funcInputs, 'I', [], @isvector);
addParameter(funcInputs, 'injectionTime', [], @isvector);

dummy.dummy = 10;
parse(funcInputs, dummy, varargin{:});
funcInputs = funcInputs.Results;
end

function checkDummyOk = checkDummy(dummy)
checkDummyOk = true;
end