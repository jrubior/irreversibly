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
    T = two_sector_alt.dynamic_resid_tt(T, y, x, params, steady_state, it_);
end
residual = zeros(35, 1);
    residual(1) = (y(7)) - (T(1)*T(2));
    residual(2) = (y(9)) - (T(3)*T(4));
    residual(3) = (y(10)) - (y(14)*T(17));
    residual(4) = (y(12)) - (y(13)+(1-params(5))*y(1));
    residual(5) = (y(14)) - (y(19)*T(5));
    residual(6) = (y(15)) - (y(7)*params(2)/y(8));
    residual(7) = (y(16)) - (y(9)*params(3)*y(15)/y(10));
    residual(8) = (y(17)) - (y(9)*(1-params(3))*y(15)/y(11));
    residual(9) = (y(18)) - (params(1)*((1-params(5))*y(44)+params(4)*y(43)*y(42)/y(12)));
    residual(10) = (log(y(19))) - (params(6)*log(y(2))+x(it_, 1));
    residual(11) = (log(y(20))) - (params(7)*log(y(3))+x(it_, 2));
    residual(12) = (y(22)) - (T(6)*T(7));
    residual(13) = (y(23)) - (y(27)*T(18));
    residual(14) = (y(25)) - (y(26)+(1-params(11))*y(4));
    residual(15) = (y(27)) - (y(32)*T(8));
    residual(16) = (y(28)) - (y(7)*params(8)/y(21));
    residual(17) = (y(29)) - (y(22)*params(9)*y(28)/y(23));
    residual(18) = (y(30)) - (y(22)*(1-params(9))*y(28)/y(24));
    residual(19) = (y(31)) - (params(1)*((1-params(11))*y(47)+params(10)*y(46)*y(45)/y(25)));
    residual(20) = (log(y(32))) - (params(12)*log(y(5))+x(it_, 3));
    residual(21) = (log(y(33))) - (params(13)*log(y(6))+x(it_, 4));
    residual(22) = (y(11)) - (T(9)*T(10));
    residual(23) = (y(13)) - (T(11)*T(12));
    residual(24) = (y(24)) - (T(13)*T(14));
    residual(25) = (y(26)) - (T(15)*T(16));
    residual(26) = (y(37)+y(36)+y(35)+y(8)+y(34)) - (y(9));
    residual(27) = (y(15)*y(34)) - (y(11)*y(17)*params(14));
    residual(28) = (y(15)*y(35)) - (y(13)*y(18)*params(15));
    residual(29) = (y(15)*y(36)) - (y(24)*y(30)*params(16));
    residual(30) = (y(15)*y(37)) - (y(26)*y(31)*params(17));
    residual(31) = (y(41)+y(40)+y(39)+y(21)+y(38)) - (y(22));
    residual(32) = (y(28)*y(38)) - (y(11)*y(17)*params(18));
    residual(33) = (y(28)*y(39)) - (y(13)*y(18)*params(19));
    residual(34) = (y(28)*y(40)) - (y(24)*y(30)*params(20));
    residual(35) = (y(28)*y(41)) - (y(26)*y(31)*params(21));

end
