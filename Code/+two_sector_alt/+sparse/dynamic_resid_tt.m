function [T_order, T] = dynamic_resid_tt(y, x, params, steady_state, T_order, T)
if T_order >= 0
    return
end
T_order = 0;
if size(T, 1) < 18
    T = [T; NaN(18 - size(T, 1), 1)];
end
T(1) = (y(37)/params(2))^params(2);
T(2) = (y(50)/params(8))^params(8);
T(3) = (y(39)/params(3))^params(3);
T(4) = (y(40)/(1-params(3)))^(1-params(3));
T(5) = (y(49)/(1-params(4)))^(1-params(4));
T(6) = (y(52)/params(9))^params(9);
T(7) = (y(53)/(1-params(9)))^(1-params(9));
T(8) = (y(62)/(1-params(10)))^(1-params(10));
T(9) = (y(63)/params(14))^params(14);
T(10) = (y(67)/params(18))^params(18);
T(11) = (y(64)/params(15))^params(15);
T(12) = (y(68)/params(19))^params(19);
T(13) = (y(65)/params(16))^params(16);
T(14) = (y(69)/params(20))^params(20);
T(15) = (y(66)/params(17))^params(17);
T(16) = (y(70)/params(21))^params(21);
T(17) = (y(6)/params(4))^params(4);
T(18) = (y(19)/params(10))^params(10);
end
