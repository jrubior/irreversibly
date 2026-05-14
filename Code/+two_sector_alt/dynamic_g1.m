function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
% function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
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
%   g1
%

if T_flag
    T = two_sector_alt.dynamic_g1_tt(T, y, x, params, steady_state, it_);
end
g1 = zeros(35, 51);
g1(1,7)=1;
g1(1,8)=(-(T(2)*1/params(2)*getPowerDeriv(y(8)/params(2),params(2),1)));
g1(1,21)=(-(T(1)*1/params(8)*getPowerDeriv(y(21)/params(8),params(8),1)));
g1(2,9)=1;
g1(2,10)=(-(T(4)*1/params(3)*getPowerDeriv(y(10)/params(3),params(3),1)));
g1(2,11)=(-(T(3)*1/(1-params(3))*getPowerDeriv(y(11)/(1-params(3)),1-params(3),1)));
g1(3,10)=1;
g1(3,1)=(-(y(14)*1/params(4)*getPowerDeriv(y(1)/params(4),params(4),1)));
g1(3,14)=(-T(17));
g1(4,1)=(-(1-params(5)));
g1(4,12)=1;
g1(4,13)=(-1);
g1(5,14)=1;
g1(5,19)=(-T(5));
g1(5,20)=(-(y(19)*1/(1-params(4))*getPowerDeriv(y(20)/(1-params(4)),1-params(4),1)));
g1(6,7)=(-(params(2)/y(8)));
g1(6,8)=(-((-(y(7)*params(2)))/(y(8)*y(8))));
g1(6,15)=1;
g1(7,9)=(-(params(3)*y(15)/y(10)));
g1(7,10)=(-((-(y(9)*params(3)*y(15)))/(y(10)*y(10))));
g1(7,15)=(-(y(9)*params(3)/y(10)));
g1(7,16)=1;
g1(8,9)=(-((1-params(3))*y(15)/y(11)));
g1(8,11)=(-((-(y(9)*(1-params(3))*y(15)))/(y(11)*y(11))));
g1(8,15)=(-(y(9)*(1-params(3))/y(11)));
g1(8,17)=1;
g1(9,42)=(-(params(1)*params(4)*y(43)/y(12)));
g1(9,12)=(-(params(1)*(-(params(4)*y(43)*y(42)))/(y(12)*y(12))));
g1(9,43)=(-(params(1)*params(4)*y(42)/y(12)));
g1(9,18)=1;
g1(9,44)=(-((1-params(5))*params(1)));
g1(10,2)=(-(params(6)*1/y(2)));
g1(10,19)=1/y(19);
g1(10,48)=(-1);
g1(11,3)=(-(params(7)*1/y(3)));
g1(11,20)=1/y(20);
g1(11,49)=(-1);
g1(12,22)=1;
g1(12,23)=(-(T(7)*1/params(9)*getPowerDeriv(y(23)/params(9),params(9),1)));
g1(12,24)=(-(T(6)*1/(1-params(9))*getPowerDeriv(y(24)/(1-params(9)),1-params(9),1)));
g1(13,23)=1;
g1(13,4)=(-(y(27)*1/params(10)*getPowerDeriv(y(4)/params(10),params(10),1)));
g1(13,27)=(-T(18));
g1(14,4)=(-(1-params(11)));
g1(14,25)=1;
g1(14,26)=(-1);
g1(15,27)=1;
g1(15,32)=(-T(8));
g1(15,33)=(-(y(32)*1/(1-params(10))*getPowerDeriv(y(33)/(1-params(10)),1-params(10),1)));
g1(16,7)=(-(params(8)/y(21)));
g1(16,21)=(-((-(y(7)*params(8)))/(y(21)*y(21))));
g1(16,28)=1;
g1(17,22)=(-(params(9)*y(28)/y(23)));
g1(17,23)=(-((-(y(22)*params(9)*y(28)))/(y(23)*y(23))));
g1(17,28)=(-(y(22)*params(9)/y(23)));
g1(17,29)=1;
g1(18,22)=(-((1-params(9))*y(28)/y(24)));
g1(18,24)=(-((-(y(22)*(1-params(9))*y(28)))/(y(24)*y(24))));
g1(18,28)=(-(y(22)*(1-params(9))/y(24)));
g1(18,30)=1;
g1(19,45)=(-(params(1)*params(10)*y(46)/y(25)));
g1(19,25)=(-(params(1)*(-(params(10)*y(46)*y(45)))/(y(25)*y(25))));
g1(19,46)=(-(params(1)*params(10)*y(45)/y(25)));
g1(19,31)=1;
g1(19,47)=(-(params(1)*(1-params(11))));
g1(20,5)=(-(params(12)*1/y(5)));
g1(20,32)=1/y(32);
g1(20,50)=(-1);
g1(21,6)=(-(params(13)*1/y(6)));
g1(21,33)=1/y(33);
g1(21,51)=(-1);
g1(22,11)=1;
g1(22,34)=(-(T(10)*1/params(14)*getPowerDeriv(y(34)/params(14),params(14),1)));
g1(22,38)=(-(T(9)*1/params(18)*getPowerDeriv(y(38)/params(18),params(18),1)));
g1(23,13)=1;
g1(23,35)=(-(T(12)*1/params(15)*getPowerDeriv(y(35)/params(15),params(15),1)));
g1(23,39)=(-(T(11)*1/params(19)*getPowerDeriv(y(39)/params(19),params(19),1)));
g1(24,24)=1;
g1(24,36)=(-(T(14)*1/params(16)*getPowerDeriv(y(36)/params(16),params(16),1)));
g1(24,40)=(-(T(13)*1/params(20)*getPowerDeriv(y(40)/params(20),params(20),1)));
g1(25,26)=1;
g1(25,37)=(-(T(16)*1/params(17)*getPowerDeriv(y(37)/params(17),params(17),1)));
g1(25,41)=(-(T(15)*1/params(21)*getPowerDeriv(y(41)/params(21),params(21),1)));
g1(26,8)=1;
g1(26,9)=(-1);
g1(26,34)=1;
g1(26,35)=1;
g1(26,36)=1;
g1(26,37)=1;
g1(27,11)=(-(y(17)*params(14)));
g1(27,15)=y(34);
g1(27,17)=(-(y(11)*params(14)));
g1(27,34)=y(15);
g1(28,13)=(-(y(18)*params(15)));
g1(28,15)=y(35);
g1(28,18)=(-(y(13)*params(15)));
g1(28,35)=y(15);
g1(29,15)=y(36);
g1(29,24)=(-(y(30)*params(16)));
g1(29,30)=(-(y(24)*params(16)));
g1(29,36)=y(15);
g1(30,15)=y(37);
g1(30,26)=(-(y(31)*params(17)));
g1(30,31)=(-(y(26)*params(17)));
g1(30,37)=y(15);
g1(31,21)=1;
g1(31,22)=(-1);
g1(31,38)=1;
g1(31,39)=1;
g1(31,40)=1;
g1(31,41)=1;
g1(32,11)=(-(y(17)*params(18)));
g1(32,17)=(-(y(11)*params(18)));
g1(32,28)=y(38);
g1(32,38)=y(28);
g1(33,13)=(-(y(18)*params(19)));
g1(33,18)=(-(y(13)*params(19)));
g1(33,28)=y(39);
g1(33,39)=y(28);
g1(34,24)=(-(y(30)*params(20)));
g1(34,28)=y(40);
g1(34,30)=(-(y(24)*params(20)));
g1(34,40)=y(28);
g1(35,26)=(-(y(31)*params(21)));
g1(35,28)=y(41);
g1(35,31)=(-(y(26)*params(21)));
g1(35,41)=y(28);

end
