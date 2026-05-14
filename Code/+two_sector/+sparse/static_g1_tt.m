function [T_order, T] = static_g1_tt(y, x, params, T_order, T)
if T_order >= 1
    return
end
[T_order, T] = two_sector.sparse.static_resid_tt(y, x, params, T_order, T);
T_order = 1;
if size(T, 1) < 26
    T = [T; NaN(26 - size(T, 1), 1)];
end
T(23) = (-(exp(y(4))*(1-params(4))*exp(y(24))/exp(y(8))));
T(24) = (-(exp(y(5))*(1-params(5))*exp(y(25))/exp(y(9))));
T(25) = (-(params(1)*exp(y(6))*params(6)*exp(y(26))/exp(y(14))));
T(26) = (-(params(1)*exp(y(7))*params(7)*exp(y(27))/exp(y(15))));
end
