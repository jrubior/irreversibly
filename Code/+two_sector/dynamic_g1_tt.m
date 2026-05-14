function T = dynamic_g1_tt(T, y, x, params, steady_state, it_)
% function T = dynamic_g1_tt(T, y, x, params, steady_state, it_)
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

assert(length(T) >= 26);

T = two_sector.dynamic_resid_tt(T, y, x, params, steady_state, it_);

T(23) = (-(exp(y(10))*(1-params(4))*exp(y(30))/exp(y(14))));
T(24) = (-(exp(y(11))*(1-params(5))*exp(y(31))/exp(y(15))));
T(25) = (-(params(1)*params(6)*exp(y(44))*exp(y(42))/exp(y(20))));
T(26) = (-(params(1)*params(7)*exp(y(45))*exp(y(43))/exp(y(21))));

end
