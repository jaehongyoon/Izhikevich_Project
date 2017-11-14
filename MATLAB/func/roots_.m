function rts = roots_(pars)
% find the fixed point of Izhikevich model for given parameters

a = pars(1); b = pars(2); c = pars(3); d = pars(4);
I = pars(5);

coeff = [.04 -5-b -140+I];
rts = roots(coeff);
end