function T = dynamic_resid_tt(T, y, x, params, steady_state, it_)
% function T = dynamic_resid_tt(T, y, x, params, steady_state, it_)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double  vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double  vector of endogenous variables in the order stored
%                                                    in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double  matrix of exogenous variables (in declaration order)
%                                                    for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double  vector of steady state values
%   params        [M_.param_nbr by 1]        double  vector of parameter values in declaration order
%   it_           scalar                     double  time period for exogenous variables for which
%                                                    to evaluate the model
%
% Output:
%   T           [#temp variables by 1]       double  vector of temporary terms
%

assert(length(T) >= 22);

T(1) = (exp(y(8))/params(2))^params(2);
T(2) = (exp(y(9))/params(3))^params(3);
T(3) = (exp(y(12))/params(4))^params(4);
T(4) = exp(y(14))/(1-params(4));
T(5) = T(4)^(1-params(4));
T(6) = (exp(y(13))/params(5))^params(5);
T(7) = exp(y(15))/(1-params(5));
T(8) = T(7)^(1-params(5));
T(9) = (exp(y(16))/params(10))^params(10);
T(10) = (exp(y(18))/params(11))^params(11);
T(11) = (exp(y(17))/params(12))^params(12);
T(12) = (exp(y(19))/params(13))^params(13);
T(13) = (exp(y(24))/params(14))^params(14);
T(14) = (exp(y(26))/params(15))^params(15);
T(15) = (exp(y(25))/params(16))^params(16);
T(16) = (exp(y(27))/params(17))^params(17);
T(17) = exp(y(40))/(1-params(6));
T(18) = exp(y(38))*T(17)^(1-params(6));
T(19) = exp(y(41))/(1-params(7));
T(20) = exp(y(39))*T(19)^(1-params(7));
T(21) = exp(y(28))*(exp(y(1))/params(6))^params(6);
T(22) = exp(y(29))*(exp(y(2))/params(7))^params(7);

end
