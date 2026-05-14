function g1 = static_g1(T, y, x, params, T_flag)
% function g1 = static_g1(T, y, x, params, T_flag)
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
%   g1
%

if T_flag
    T = two_sector_alt.static_g1_tt(T, y, x, params);
end
g1 = zeros(35, 35);
g1(1,1)=1;
g1(1,2)=(-(T(2)*1/params(2)*getPowerDeriv(y(2)/params(2),params(2),1)));
g1(1,15)=(-(T(1)*1/params(8)*getPowerDeriv(y(15)/params(8),params(8),1)));
g1(2,3)=1;
g1(2,4)=(-(T(4)*1/params(3)*getPowerDeriv(y(4)/params(3),params(3),1)));
g1(2,5)=(-(T(3)*1/(1-params(3))*getPowerDeriv(y(5)/(1-params(3)),1-params(3),1)));
g1(3,4)=1;
g1(3,6)=(-(y(8)*1/params(4)*getPowerDeriv(y(6)/params(4),params(4),1)));
g1(3,8)=(-T(5));
g1(4,6)=1-(1-params(5));
g1(4,7)=(-1);
g1(5,8)=1;
g1(5,13)=(-T(6));
g1(5,14)=(-(y(13)*1/(1-params(4))*getPowerDeriv(y(14)/(1-params(4)),1-params(4),1)));
g1(6,1)=(-(params(2)/y(2)));
g1(6,2)=(-((-(y(1)*params(2)))/(y(2)*y(2))));
g1(6,9)=1;
g1(7,3)=(-(params(3)*y(9)/y(4)));
g1(7,4)=(-((-(y(3)*params(3)*y(9)))/(y(4)*y(4))));
g1(7,9)=(-(y(3)*params(3)/y(4)));
g1(7,10)=1;
g1(8,3)=(-((1-params(3))*y(9)/y(5)));
g1(8,5)=(-((-(y(3)*(1-params(3))*y(9)))/(y(5)*y(5))));
g1(8,9)=(-(y(3)*(1-params(3))/y(5)));
g1(8,11)=1;
g1(9,4)=(-(params(1)*params(4)*y(10)/y(6)));
g1(9,6)=(-(params(1)*(-(y(4)*params(4)*y(10)))/(y(6)*y(6))));
g1(9,10)=(-(params(1)*y(4)*params(4)/y(6)));
g1(9,12)=1-(1-params(5))*params(1);
g1(10,13)=1/y(13)-params(6)*1/y(13);
g1(11,14)=1/y(14)-params(7)*1/y(14);
g1(12,16)=1;
g1(12,17)=(-(T(8)*1/params(9)*getPowerDeriv(y(17)/params(9),params(9),1)));
g1(12,18)=(-(T(7)*1/(1-params(9))*getPowerDeriv(y(18)/(1-params(9)),1-params(9),1)));
g1(13,17)=1;
g1(13,19)=(-(y(21)*1/params(10)*getPowerDeriv(y(19)/params(10),params(10),1)));
g1(13,21)=(-T(9));
g1(14,19)=1-(1-params(11));
g1(14,20)=(-1);
g1(15,21)=1;
g1(15,26)=(-T(10));
g1(15,27)=(-(y(26)*1/(1-params(10))*getPowerDeriv(y(27)/(1-params(10)),1-params(10),1)));
g1(16,1)=(-(params(8)/y(15)));
g1(16,15)=(-((-(y(1)*params(8)))/(y(15)*y(15))));
g1(16,22)=1;
g1(17,16)=(-(params(9)*y(22)/y(17)));
g1(17,17)=(-((-(y(16)*params(9)*y(22)))/(y(17)*y(17))));
g1(17,22)=(-(y(16)*params(9)/y(17)));
g1(17,23)=1;
g1(18,16)=(-((1-params(9))*y(22)/y(18)));
g1(18,18)=(-((-(y(16)*(1-params(9))*y(22)))/(y(18)*y(18))));
g1(18,22)=(-(y(16)*(1-params(9))/y(18)));
g1(18,24)=1;
g1(19,17)=(-(params(1)*params(10)*y(23)/y(19)));
g1(19,19)=(-(params(1)*(-(y(17)*params(10)*y(23)))/(y(19)*y(19))));
g1(19,23)=(-(params(1)*y(17)*params(10)/y(19)));
g1(19,25)=1-params(1)*(1-params(11));
g1(20,26)=1/y(26)-params(12)*1/y(26);
g1(21,27)=1/y(27)-params(13)*1/y(27);
g1(22,5)=1;
g1(22,28)=(-(T(12)*1/params(14)*getPowerDeriv(y(28)/params(14),params(14),1)));
g1(22,32)=(-(T(11)*1/params(18)*getPowerDeriv(y(32)/params(18),params(18),1)));
g1(23,7)=1;
g1(23,29)=(-(T(14)*1/params(15)*getPowerDeriv(y(29)/params(15),params(15),1)));
g1(23,33)=(-(T(13)*1/params(19)*getPowerDeriv(y(33)/params(19),params(19),1)));
g1(24,18)=1;
g1(24,30)=(-(T(16)*1/params(16)*getPowerDeriv(y(30)/params(16),params(16),1)));
g1(24,34)=(-(T(15)*1/params(20)*getPowerDeriv(y(34)/params(20),params(20),1)));
g1(25,20)=1;
g1(25,31)=(-(T(18)*1/params(17)*getPowerDeriv(y(31)/params(17),params(17),1)));
g1(25,35)=(-(T(17)*1/params(21)*getPowerDeriv(y(35)/params(21),params(21),1)));
g1(26,2)=1;
g1(26,3)=(-1);
g1(26,28)=1;
g1(26,29)=1;
g1(26,30)=1;
g1(26,31)=1;
g1(27,5)=(-(y(11)*params(14)));
g1(27,9)=y(28);
g1(27,11)=(-(y(5)*params(14)));
g1(27,28)=y(9);
g1(28,7)=(-(y(12)*params(15)));
g1(28,9)=y(29);
g1(28,12)=(-(y(7)*params(15)));
g1(28,29)=y(9);
g1(29,9)=y(30);
g1(29,18)=(-(y(24)*params(16)));
g1(29,24)=(-(y(18)*params(16)));
g1(29,30)=y(9);
g1(30,9)=y(31);
g1(30,20)=(-(y(25)*params(17)));
g1(30,25)=(-(y(20)*params(17)));
g1(30,31)=y(9);
g1(31,15)=1;
g1(31,16)=(-1);
g1(31,32)=1;
g1(31,33)=1;
g1(31,34)=1;
g1(31,35)=1;
g1(32,5)=(-(y(11)*params(18)));
g1(32,11)=(-(y(5)*params(18)));
g1(32,22)=y(32);
g1(32,32)=y(22);
g1(33,7)=(-(y(12)*params(19)));
g1(33,12)=(-(y(7)*params(19)));
g1(33,22)=y(33);
g1(33,33)=y(22);
g1(34,18)=(-(y(24)*params(20)));
g1(34,22)=y(34);
g1(34,24)=(-(y(18)*params(20)));
g1(34,34)=y(22);
g1(35,20)=(-(y(25)*params(21)));
g1(35,22)=y(35);
g1(35,25)=(-(y(20)*params(21)));
g1(35,35)=y(22);

end
