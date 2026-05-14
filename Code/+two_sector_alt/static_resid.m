function residual = static_resid(T, y, x, params, T_flag)
% function residual = static_resid(T, y, x, params, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%                                              to evaluate the model
%   T_flag    boolean                 boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   residual
%

if T_flag
    T = two_sector_alt.static_resid_tt(T, y, x, params);
end
residual = zeros(35, 1);
    residual(1) = (y(1)) - (T(1)*T(2));
    residual(2) = (y(3)) - (T(3)*T(4));
    residual(3) = (y(4)) - (y(8)*T(5));
    residual(4) = (y(6)) - (y(7)+y(6)*(1-params(5)));
    residual(5) = (y(8)) - (y(13)*T(6));
    residual(6) = (y(9)) - (y(1)*params(2)/y(2));
    residual(7) = (y(10)) - (y(3)*params(3)*y(9)/y(4));
    residual(8) = (y(11)) - (y(3)*(1-params(3))*y(9)/y(5));
    residual(9) = (y(12)) - (params(1)*((1-params(5))*y(12)+y(4)*params(4)*y(10)/y(6)));
    residual(10) = (log(y(13))) - (log(y(13))*params(6)+x(1));
    residual(11) = (log(y(14))) - (log(y(14))*params(7)+x(2));
    residual(12) = (y(16)) - (T(7)*T(8));
    residual(13) = (y(17)) - (y(21)*T(9));
    residual(14) = (y(19)) - (y(20)+y(19)*(1-params(11)));
    residual(15) = (y(21)) - (y(26)*T(10));
    residual(16) = (y(22)) - (y(1)*params(8)/y(15));
    residual(17) = (y(23)) - (y(16)*params(9)*y(22)/y(17));
    residual(18) = (y(24)) - (y(16)*(1-params(9))*y(22)/y(18));
    residual(19) = (y(25)) - (params(1)*((1-params(11))*y(25)+y(17)*params(10)*y(23)/y(19)));
    residual(20) = (log(y(26))) - (log(y(26))*params(12)+x(3));
    residual(21) = (log(y(27))) - (log(y(27))*params(13)+x(4));
    residual(22) = (y(5)) - (T(11)*T(12));
    residual(23) = (y(7)) - (T(13)*T(14));
    residual(24) = (y(18)) - (T(15)*T(16));
    residual(25) = (y(20)) - (T(17)*T(18));
    residual(26) = (y(31)+y(30)+y(29)+y(2)+y(28)) - (y(3));
    residual(27) = (y(9)*y(28)) - (y(5)*y(11)*params(14));
    residual(28) = (y(9)*y(29)) - (y(7)*y(12)*params(15));
    residual(29) = (y(9)*y(30)) - (y(18)*y(24)*params(16));
    residual(30) = (y(9)*y(31)) - (y(20)*y(25)*params(17));
    residual(31) = (y(35)+y(34)+y(33)+y(15)+y(32)) - (y(16));
    residual(32) = (y(22)*y(32)) - (y(5)*y(11)*params(18));
    residual(33) = (y(22)*y(33)) - (y(7)*y(12)*params(19));
    residual(34) = (y(22)*y(34)) - (y(18)*y(24)*params(20));
    residual(35) = (y(22)*y(35)) - (y(20)*y(25)*params(21));

end
