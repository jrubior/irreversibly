function [T_order, T] = static_resid_tt(y, x, params, T_order, T)
if T_order >= 0
    return
end
T_order = 0;
if size(T, 1) < 22
    T = [T; NaN(22 - size(T, 1), 1)];
end
T(1) = (exp(y(2))/params(2))^params(2);
T(2) = (exp(y(3))/params(3))^params(3);
T(3) = (exp(y(6))/params(4))^params(4);
T(4) = exp(y(8))/(1-params(4));
T(5) = T(4)^(1-params(4));
T(6) = (exp(y(7))/params(5))^params(5);
T(7) = exp(y(9))/(1-params(5));
T(8) = T(7)^(1-params(5));
T(9) = (exp(y(10))/params(10))^params(10);
T(10) = (exp(y(12))/params(11))^params(11);
T(11) = (exp(y(11))/params(12))^params(12);
T(12) = (exp(y(13))/params(13))^params(13);
T(13) = exp(y(22))*(exp(y(14))/params(6))^params(6);
T(14) = exp(y(23))*(exp(y(15))/params(7))^params(7);
T(15) = (exp(y(18))/params(14))^params(14);
T(16) = (exp(y(20))/params(15))^params(15);
T(17) = (exp(y(19))/params(16))^params(16);
T(18) = (exp(y(21))/params(17))^params(17);
T(19) = exp(y(34))/(1-params(6));
T(20) = exp(y(32))*T(19)^(1-params(6));
T(21) = exp(y(35))/(1-params(7));
T(22) = exp(y(33))*T(21)^(1-params(7));
end
