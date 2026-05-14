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

assert(length(T) >= 18);

T(1) = (y(8)/params(2))^params(2);
T(2) = (y(21)/params(8))^params(8);
T(3) = (y(10)/params(3))^params(3);
T(4) = (y(11)/(1-params(3)))^(1-params(3));
T(5) = (y(20)/(1-params(4)))^(1-params(4));
T(6) = (y(23)/params(9))^params(9);
T(7) = (y(24)/(1-params(9)))^(1-params(9));
T(8) = (y(33)/(1-params(10)))^(1-params(10));
T(9) = (y(34)/params(14))^params(14);
T(10) = (y(38)/params(18))^params(18);
T(11) = (y(35)/params(15))^params(15);
T(12) = (y(39)/params(19))^params(19);
T(13) = (y(36)/params(16))^params(16);
T(14) = (y(40)/params(20))^params(20);
T(15) = (y(37)/params(17))^params(17);
T(16) = (y(41)/params(21))^params(21);
T(17) = (y(1)/params(4))^params(4);
T(18) = (y(4)/params(10))^params(10);

end
