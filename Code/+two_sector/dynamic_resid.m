function residual = dynamic_resid(T, y, x, params, steady_state, it_, T_flag)
% function residual = dynamic_resid(T, y, x, params, steady_state, it_, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double   vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double   vector of endogenous variables in the order stored
%                                                     in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double   matrix of exogenous variables (in declaration order)
%                                                     for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double   vector of steady state values
%   params        [M_.param_nbr by 1]        double   vector of parameter values in declaration order
%   it_           scalar                     double   time period for exogenous variables for which
%                                                     to evaluate the model
%   T_flag        boolean                    boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   residual
%

if T_flag
    T = two_sector.dynamic_resid_tt(T, y, x, params, steady_state, it_);
end
residual = zeros(35, 1);
    residual(1) = (exp(y(7))) - (T(1)*T(2));
    residual(2) = (exp(y(10))) - (T(3)*T(5));
    residual(3) = (exp(y(11))) - (T(6)*T(8));
    residual(4) = (exp(y(14))) - (T(9)*T(10));
    residual(5) = (exp(y(15))) - (T(11)*T(12));
    residual(6) = (exp(y(12))) - (T(21));
    residual(7) = (exp(y(13))) - (T(22));
    residual(8) = (exp(y(20))) - (exp(y(22))+(1-params(8))*exp(y(1)));
    residual(9) = (exp(y(21))) - (exp(y(23))+(1-params(9))*exp(y(2)));
    residual(10) = (exp(y(22))) - (T(13)*T(14));
    residual(11) = (exp(y(23))) - (T(15)*T(16));
    residual(12) = (exp(y(25))+exp(y(24))+exp(y(17))+exp(y(8))+exp(y(16))) - (exp(y(10)));
    residual(13) = (exp(y(27))+exp(y(26))+exp(y(19))+exp(y(9))+exp(y(18))) - (exp(y(11)));
    residual(14) = (exp(y(28))) - (T(18));
    residual(15) = (exp(y(29))) - (T(20));
    residual(16) = (exp(y(30))) - (exp(y(7))*params(2)/exp(y(8)));
    residual(17) = (exp(y(31))) - (exp(y(7))*params(3)/exp(y(9)));
    residual(18) = (exp(y(32))) - (exp(y(10))*params(4)*exp(y(30))/exp(y(12)));
    residual(19) = (exp(y(33))) - (exp(y(11))*params(5)*exp(y(31))/exp(y(13)));
    residual(20) = (exp(y(34))) - (exp(y(10))*(1-params(4))*exp(y(30))/exp(y(14)));
    residual(21) = (exp(y(35))) - (exp(y(11))*(1-params(5))*exp(y(31))/exp(y(15)));
    residual(22) = (exp(y(16))*exp(y(30))) - (exp(y(14))*params(10)*exp(y(34)));
    residual(23) = (exp(y(18))*exp(y(31))) - (exp(y(14))*params(11)*exp(y(34)));
    residual(24) = (exp(y(17))*exp(y(30))) - (exp(y(15))*params(12)*exp(y(35)));
    residual(25) = (exp(y(19))*exp(y(31))) - (exp(y(15))*params(13)*exp(y(35)));
    residual(26) = (exp(y(24))*exp(y(30))) - (exp(y(22))*params(14)*exp(y(36)));
    residual(27) = (exp(y(26))*exp(y(31))) - (exp(y(22))*params(15)*exp(y(36)));
    residual(28) = (exp(y(25))*exp(y(30))) - (exp(y(23))*params(16)*exp(y(37)));
    residual(29) = (exp(y(27))*exp(y(31))) - (exp(y(23))*params(17)*exp(y(37)));
    residual(30) = (exp(y(36))) - (params(1)*((1-params(8))*exp(y(46))+params(6)*exp(y(44))*exp(y(42))/exp(y(20))));
    residual(31) = (exp(y(37))) - (params(1)*((1-params(9))*exp(y(47))+params(7)*exp(y(45))*exp(y(43))/exp(y(21))));
    residual(32) = (y(38)) - (params(18)*y(3)+x(it_, 1));
    residual(33) = (y(39)) - (params(19)*y(4)+x(it_, 2));
    residual(34) = (y(40)) - (params(20)*y(5)+x(it_, 3));
    residual(35) = (y(41)) - (params(21)*y(6)+x(it_, 4));

end
