function T = static_resid_tt(T, y, x, params)
% function T = static_resid_tt(T, y, x, params)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%
% Output:
%   T         [#temp variables by 1]  double   vector of temporary terms
%

assert(length(T) >= 18);

T(1) = (y(2)/params(2))^params(2);
T(2) = (y(15)/params(8))^params(8);
T(3) = (y(4)/params(3))^params(3);
T(4) = (y(5)/(1-params(3)))^(1-params(3));
T(5) = (y(6)/params(4))^params(4);
T(6) = (y(14)/(1-params(4)))^(1-params(4));
T(7) = (y(17)/params(9))^params(9);
T(8) = (y(18)/(1-params(9)))^(1-params(9));
T(9) = (y(19)/params(10))^params(10);
T(10) = (y(27)/(1-params(10)))^(1-params(10));
T(11) = (y(28)/params(14))^params(14);
T(12) = (y(32)/params(18))^params(18);
T(13) = (y(29)/params(15))^params(15);
T(14) = (y(33)/params(19))^params(19);
T(15) = (y(30)/params(16))^params(16);
T(16) = (y(34)/params(20))^params(20);
T(17) = (y(31)/params(17))^params(17);
T(18) = (y(35)/params(21))^params(21);

end
