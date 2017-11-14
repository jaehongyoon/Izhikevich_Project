function J = Izhikevich_Jacobian(v, u, pars)
% Jacobian matrix of Izhikevich model
fv =@(v,u) .08*v - 5;
fu =@(v,u) -1;
gv =@(v,u) pars(1)*pars(2);
gu =@(v,u) -pars(1);

J = [fv(v,u) fu(v,u);gv(v,u) gu(v,u)];
end