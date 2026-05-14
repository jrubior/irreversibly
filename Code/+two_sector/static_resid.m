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
    T = two_sector.static_resid_tt(T, y, x, params);
end
residual = zeros(35, 1);
    residual(1) = (exp(y(1))) - (T(1)*T(2));
    residual(2) = (exp(y(4))) - (T(3)*T(5));
    residual(3) = (exp(y(5))) - (T(6)*T(8));
    residual(4) = (exp(y(8))) - (T(9)*T(10));
    residual(5) = (exp(y(9))) - (T(11)*T(12));
    residual(6) = (exp(y(6))) - (T(13));
    residual(7) = (exp(y(7))) - (T(14));
    residual(8) = (exp(y(14))) - (exp(y(16))+exp(y(14))*(1-params(8)));
    residual(9) = (exp(y(15))) - (exp(y(17))+exp(y(15))*(1-params(9)));
    residual(10) = (exp(y(16))) - (T(15)*T(16));
    residual(11) = (exp(y(17))) - (T(17)*T(18));
    residual(12) = (exp(y(19))+exp(y(18))+exp(y(11))+exp(y(2))+exp(y(10))) - (exp(y(4)));
    residual(13) = (exp(y(21))+exp(y(20))+exp(y(13))+exp(y(3))+exp(y(12))) - (exp(y(5)));
    residual(14) = (exp(y(22))) - (T(20));
    residual(15) = (exp(y(23))) - (T(22));
    residual(16) = (exp(y(24))) - (exp(y(1))*params(2)/exp(y(2)));
    residual(17) = (exp(y(25))) - (exp(y(1))*params(3)/exp(y(3)));
    residual(18) = (exp(y(26))) - (exp(y(4))*params(4)*exp(y(24))/exp(y(6)));
    residual(19) = (exp(y(27))) - (exp(y(5))*params(5)*exp(y(25))/exp(y(7)));
    residual(20) = (exp(y(28))) - (exp(y(4))*(1-params(4))*exp(y(24))/exp(y(8)));
    residual(21) = (exp(y(29))) - (exp(y(5))*(1-params(5))*exp(y(25))/exp(y(9)));
    residual(22) = (exp(y(10))*exp(y(24))) - (exp(y(8))*params(10)*exp(y(28)));
    residual(23) = (exp(y(12))*exp(y(25))) - (exp(y(8))*params(11)*exp(y(28)));
    residual(24) = (exp(y(11))*exp(y(24))) - (exp(y(9))*params(12)*exp(y(29)));
    residual(25) = (exp(y(13))*exp(y(25))) - (exp(y(9))*params(13)*exp(y(29)));
    residual(26) = (exp(y(18))*exp(y(24))) - (exp(y(16))*params(14)*exp(y(30)));
    residual(27) = (exp(y(20))*exp(y(25))) - (exp(y(16))*params(15)*exp(y(30)));
    residual(28) = (exp(y(19))*exp(y(24))) - (exp(y(17))*params(16)*exp(y(31)));
    residual(29) = (exp(y(21))*exp(y(25))) - (exp(y(17))*params(17)*exp(y(31)));
    residual(30) = (exp(y(30))) - (params(1)*((1-params(8))*exp(y(30))+exp(y(6))*params(6)*exp(y(26))/exp(y(14))));
    residual(31) = (exp(y(31))) - (params(1)*((1-params(9))*exp(y(31))+exp(y(7))*params(7)*exp(y(27))/exp(y(15))));
    residual(32) = (y(32)) - (y(32)*params(18)+x(1));
    residual(33) = (y(33)) - (y(33)*params(19)+x(2));
    residual(34) = (y(34)) - (y(34)*params(20)+x(3));
    residual(35) = (y(35)) - (y(35)*params(21)+x(4));

end
