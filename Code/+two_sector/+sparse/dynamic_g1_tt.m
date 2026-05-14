function [T_order, T] = dynamic_g1_tt(y, x, params, steady_state, T_order, T)
if T_order >= 1
    return
end
[T_order, T] = two_sector.sparse.dynamic_resid_tt(y, x, params, steady_state, T_order, T);
T_order = 1;
if size(T, 1) < 26
    T = [T; NaN(26 - size(T, 1), 1)];
end
T(23) = (-(exp(y(39))*(1-params(4))*exp(y(59))/exp(y(43))));
T(24) = (-(exp(y(40))*(1-params(5))*exp(y(60))/exp(y(44))));
T(25) = (-(params(1)*params(6)*exp(y(96))*exp(y(76))/exp(y(49))));
T(26) = (-(params(1)*params(7)*exp(y(97))*exp(y(77))/exp(y(50))));
end
