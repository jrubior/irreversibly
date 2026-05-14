function T = static_g1_tt(T, y, x, params)
% function T = static_g1_tt(T, y, x, params)
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

assert(length(T) >= 26);

T = two_sector.static_resid_tt(T, y, x, params);

T(23) = (-(exp(y(4))*(1-params(4))*exp(y(24))/exp(y(8))));
T(24) = (-(exp(y(5))*(1-params(5))*exp(y(25))/exp(y(9))));
T(25) = (-(params(1)*exp(y(6))*params(6)*exp(y(26))/exp(y(14))));
T(26) = (-(params(1)*exp(y(7))*params(7)*exp(y(27))/exp(y(15))));

end
