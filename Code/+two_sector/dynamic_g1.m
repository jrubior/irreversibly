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
    T = two_sector.dynamic_g1_tt(T, y, x, params, steady_state, it_);
end
g1 = zeros(35, 51);
g1(1,7)=exp(y(7));
g1(1,8)=(-(T(2)*exp(y(8))/params(2)*getPowerDeriv(exp(y(8))/params(2),params(2),1)));
g1(1,9)=(-(T(1)*exp(y(9))/params(3)*getPowerDeriv(exp(y(9))/params(3),params(3),1)));
g1(2,10)=exp(y(10));
g1(2,12)=(-(T(5)*exp(y(12))/params(4)*getPowerDeriv(exp(y(12))/params(4),params(4),1)));
g1(2,14)=(-(T(3)*T(4)*getPowerDeriv(T(4),1-params(4),1)));
g1(3,11)=exp(y(11));
g1(3,13)=(-(T(8)*exp(y(13))/params(5)*getPowerDeriv(exp(y(13))/params(5),params(5),1)));
g1(3,15)=(-(T(6)*T(7)*getPowerDeriv(T(7),1-params(5),1)));
g1(4,14)=exp(y(14));
g1(4,16)=(-(T(10)*exp(y(16))/params(10)*getPowerDeriv(exp(y(16))/params(10),params(10),1)));
g1(4,18)=(-(T(9)*exp(y(18))/params(11)*getPowerDeriv(exp(y(18))/params(11),params(11),1)));
g1(5,15)=exp(y(15));
g1(5,17)=(-(T(12)*exp(y(17))/params(12)*getPowerDeriv(exp(y(17))/params(12),params(12),1)));
g1(5,19)=(-(T(11)*exp(y(19))/params(13)*getPowerDeriv(exp(y(19))/params(13),params(13),1)));
g1(6,12)=exp(y(12));
g1(6,1)=(-(exp(y(28))*exp(y(1))/params(6)*getPowerDeriv(exp(y(1))/params(6),params(6),1)));
g1(6,28)=(-T(21));
g1(7,13)=exp(y(13));
g1(7,2)=(-(exp(y(29))*exp(y(2))/params(7)*getPowerDeriv(exp(y(2))/params(7),params(7),1)));
g1(7,29)=(-T(22));
g1(8,1)=(-((1-params(8))*exp(y(1))));
g1(8,20)=exp(y(20));
g1(8,22)=(-exp(y(22)));
g1(9,2)=(-((1-params(9))*exp(y(2))));
g1(9,21)=exp(y(21));
g1(9,23)=(-exp(y(23)));
g1(10,22)=exp(y(22));
g1(10,24)=(-(T(14)*exp(y(24))/params(14)*getPowerDeriv(exp(y(24))/params(14),params(14),1)));
g1(10,26)=(-(T(13)*exp(y(26))/params(15)*getPowerDeriv(exp(y(26))/params(15),params(15),1)));
g1(11,23)=exp(y(23));
g1(11,25)=(-(T(16)*exp(y(25))/params(16)*getPowerDeriv(exp(y(25))/params(16),params(16),1)));
g1(11,27)=(-(T(15)*exp(y(27))/params(17)*getPowerDeriv(exp(y(27))/params(17),params(17),1)));
g1(12,8)=exp(y(8));
g1(12,10)=(-exp(y(10)));
g1(12,16)=exp(y(16));
g1(12,17)=exp(y(17));
g1(12,24)=exp(y(24));
g1(12,25)=exp(y(25));
g1(13,9)=exp(y(9));
g1(13,11)=(-exp(y(11)));
g1(13,18)=exp(y(18));
g1(13,19)=exp(y(19));
g1(13,26)=exp(y(26));
g1(13,27)=exp(y(27));
g1(14,28)=exp(y(28));
g1(14,38)=(-T(18));
g1(14,40)=(-(exp(y(38))*T(17)*getPowerDeriv(T(17),1-params(6),1)));
g1(15,29)=exp(y(29));
g1(15,39)=(-T(20));
g1(15,41)=(-(exp(y(39))*T(19)*getPowerDeriv(T(19),1-params(7),1)));
g1(16,7)=(-(exp(y(7))*params(2)/exp(y(8))));
g1(16,8)=(-((-(exp(y(8))*exp(y(7))*params(2)))/(exp(y(8))*exp(y(8)))));
g1(16,30)=exp(y(30));
g1(17,7)=(-(exp(y(7))*params(3)/exp(y(9))));
g1(17,9)=(-((-(exp(y(9))*exp(y(7))*params(3)))/(exp(y(9))*exp(y(9)))));
g1(17,31)=exp(y(31));
g1(18,10)=(-(exp(y(10))*params(4)*exp(y(30))/exp(y(12))));
g1(18,12)=(-((-(exp(y(12))*exp(y(10))*params(4)*exp(y(30))))/(exp(y(12))*exp(y(12)))));
g1(18,30)=(-(exp(y(10))*params(4)*exp(y(30))/exp(y(12))));
g1(18,32)=exp(y(32));
g1(19,11)=(-(exp(y(11))*params(5)*exp(y(31))/exp(y(13))));
g1(19,13)=(-((-(exp(y(13))*exp(y(11))*params(5)*exp(y(31))))/(exp(y(13))*exp(y(13)))));
g1(19,31)=(-(exp(y(11))*params(5)*exp(y(31))/exp(y(13))));
g1(19,33)=exp(y(33));
g1(20,10)=T(23);
g1(20,14)=(-((-(exp(y(14))*exp(y(10))*(1-params(4))*exp(y(30))))/(exp(y(14))*exp(y(14)))));
g1(20,30)=T(23);
g1(20,34)=exp(y(34));
g1(21,11)=T(24);
g1(21,15)=(-((-(exp(y(15))*exp(y(11))*(1-params(5))*exp(y(31))))/(exp(y(15))*exp(y(15)))));
g1(21,31)=T(24);
g1(21,35)=exp(y(35));
g1(22,14)=(-(exp(y(14))*params(10)*exp(y(34))));
g1(22,16)=exp(y(16))*exp(y(30));
g1(22,30)=exp(y(16))*exp(y(30));
g1(22,34)=(-(exp(y(14))*params(10)*exp(y(34))));
g1(23,14)=(-(exp(y(14))*params(11)*exp(y(34))));
g1(23,18)=exp(y(18))*exp(y(31));
g1(23,31)=exp(y(18))*exp(y(31));
g1(23,34)=(-(exp(y(14))*params(11)*exp(y(34))));
g1(24,15)=(-(exp(y(15))*params(12)*exp(y(35))));
g1(24,17)=exp(y(17))*exp(y(30));
g1(24,30)=exp(y(17))*exp(y(30));
g1(24,35)=(-(exp(y(15))*params(12)*exp(y(35))));
g1(25,15)=(-(exp(y(15))*params(13)*exp(y(35))));
g1(25,19)=exp(y(19))*exp(y(31));
g1(25,31)=exp(y(19))*exp(y(31));
g1(25,35)=(-(exp(y(15))*params(13)*exp(y(35))));
g1(26,22)=(-(exp(y(22))*params(14)*exp(y(36))));
g1(26,24)=exp(y(24))*exp(y(30));
g1(26,30)=exp(y(24))*exp(y(30));
g1(26,36)=(-(exp(y(22))*params(14)*exp(y(36))));
g1(27,22)=(-(exp(y(22))*params(15)*exp(y(36))));
g1(27,26)=exp(y(26))*exp(y(31));
g1(27,31)=exp(y(26))*exp(y(31));
g1(27,36)=(-(exp(y(22))*params(15)*exp(y(36))));
g1(28,23)=(-(exp(y(23))*params(16)*exp(y(37))));
g1(28,25)=exp(y(25))*exp(y(30));
g1(28,30)=exp(y(25))*exp(y(30));
g1(28,37)=(-(exp(y(23))*params(16)*exp(y(37))));
g1(29,23)=(-(exp(y(23))*params(17)*exp(y(37))));
g1(29,27)=exp(y(27))*exp(y(31));
g1(29,31)=exp(y(27))*exp(y(31));
g1(29,37)=(-(exp(y(23))*params(17)*exp(y(37))));
g1(30,42)=T(25);
g1(30,20)=(-(params(1)*(-(exp(y(20))*params(6)*exp(y(44))*exp(y(42))))/(exp(y(20))*exp(y(20)))));
g1(30,44)=T(25);
g1(30,36)=exp(y(36));
g1(30,46)=(-(params(1)*(1-params(8))*exp(y(46))));
g1(31,43)=T(26);
g1(31,21)=(-(params(1)*(-(exp(y(21))*params(7)*exp(y(45))*exp(y(43))))/(exp(y(21))*exp(y(21)))));
g1(31,45)=T(26);
g1(31,37)=exp(y(37));
g1(31,47)=(-(params(1)*(1-params(9))*exp(y(47))));
g1(32,3)=(-params(18));
g1(32,38)=1;
g1(32,48)=(-1);
g1(33,4)=(-params(19));
g1(33,39)=1;
g1(33,49)=(-1);
g1(34,5)=(-params(20));
g1(34,40)=1;
g1(34,50)=(-1);
g1(35,6)=(-params(21));
g1(35,41)=1;
g1(35,51)=(-1);

end
