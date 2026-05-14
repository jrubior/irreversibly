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
    T = two_sector.static_g1_tt(T, y, x, params);
end
g1 = zeros(35, 35);
g1(1,1)=exp(y(1));
g1(1,2)=(-(T(2)*exp(y(2))/params(2)*getPowerDeriv(exp(y(2))/params(2),params(2),1)));
g1(1,3)=(-(T(1)*exp(y(3))/params(3)*getPowerDeriv(exp(y(3))/params(3),params(3),1)));
g1(2,4)=exp(y(4));
g1(2,6)=(-(T(5)*exp(y(6))/params(4)*getPowerDeriv(exp(y(6))/params(4),params(4),1)));
g1(2,8)=(-(T(3)*T(4)*getPowerDeriv(T(4),1-params(4),1)));
g1(3,5)=exp(y(5));
g1(3,7)=(-(T(8)*exp(y(7))/params(5)*getPowerDeriv(exp(y(7))/params(5),params(5),1)));
g1(3,9)=(-(T(6)*T(7)*getPowerDeriv(T(7),1-params(5),1)));
g1(4,8)=exp(y(8));
g1(4,10)=(-(T(10)*exp(y(10))/params(10)*getPowerDeriv(exp(y(10))/params(10),params(10),1)));
g1(4,12)=(-(T(9)*exp(y(12))/params(11)*getPowerDeriv(exp(y(12))/params(11),params(11),1)));
g1(5,9)=exp(y(9));
g1(5,11)=(-(T(12)*exp(y(11))/params(12)*getPowerDeriv(exp(y(11))/params(12),params(12),1)));
g1(5,13)=(-(T(11)*exp(y(13))/params(13)*getPowerDeriv(exp(y(13))/params(13),params(13),1)));
g1(6,6)=exp(y(6));
g1(6,14)=(-(exp(y(22))*exp(y(14))/params(6)*getPowerDeriv(exp(y(14))/params(6),params(6),1)));
g1(6,22)=(-T(13));
g1(7,7)=exp(y(7));
g1(7,15)=(-(exp(y(23))*exp(y(15))/params(7)*getPowerDeriv(exp(y(15))/params(7),params(7),1)));
g1(7,23)=(-T(14));
g1(8,14)=exp(y(14))-exp(y(14))*(1-params(8));
g1(8,16)=(-exp(y(16)));
g1(9,15)=exp(y(15))-exp(y(15))*(1-params(9));
g1(9,17)=(-exp(y(17)));
g1(10,16)=exp(y(16));
g1(10,18)=(-(T(16)*exp(y(18))/params(14)*getPowerDeriv(exp(y(18))/params(14),params(14),1)));
g1(10,20)=(-(T(15)*exp(y(20))/params(15)*getPowerDeriv(exp(y(20))/params(15),params(15),1)));
g1(11,17)=exp(y(17));
g1(11,19)=(-(T(18)*exp(y(19))/params(16)*getPowerDeriv(exp(y(19))/params(16),params(16),1)));
g1(11,21)=(-(T(17)*exp(y(21))/params(17)*getPowerDeriv(exp(y(21))/params(17),params(17),1)));
g1(12,2)=exp(y(2));
g1(12,4)=(-exp(y(4)));
g1(12,10)=exp(y(10));
g1(12,11)=exp(y(11));
g1(12,18)=exp(y(18));
g1(12,19)=exp(y(19));
g1(13,3)=exp(y(3));
g1(13,5)=(-exp(y(5)));
g1(13,12)=exp(y(12));
g1(13,13)=exp(y(13));
g1(13,20)=exp(y(20));
g1(13,21)=exp(y(21));
g1(14,22)=exp(y(22));
g1(14,32)=(-T(20));
g1(14,34)=(-(exp(y(32))*T(19)*getPowerDeriv(T(19),1-params(6),1)));
g1(15,23)=exp(y(23));
g1(15,33)=(-T(22));
g1(15,35)=(-(exp(y(33))*T(21)*getPowerDeriv(T(21),1-params(7),1)));
g1(16,1)=(-(exp(y(1))*params(2)/exp(y(2))));
g1(16,2)=(-((-(exp(y(2))*exp(y(1))*params(2)))/(exp(y(2))*exp(y(2)))));
g1(16,24)=exp(y(24));
g1(17,1)=(-(exp(y(1))*params(3)/exp(y(3))));
g1(17,3)=(-((-(exp(y(3))*exp(y(1))*params(3)))/(exp(y(3))*exp(y(3)))));
g1(17,25)=exp(y(25));
g1(18,4)=(-(exp(y(4))*params(4)*exp(y(24))/exp(y(6))));
g1(18,6)=(-((-(exp(y(6))*exp(y(4))*params(4)*exp(y(24))))/(exp(y(6))*exp(y(6)))));
g1(18,24)=(-(exp(y(4))*params(4)*exp(y(24))/exp(y(6))));
g1(18,26)=exp(y(26));
g1(19,5)=(-(exp(y(5))*params(5)*exp(y(25))/exp(y(7))));
g1(19,7)=(-((-(exp(y(7))*exp(y(5))*params(5)*exp(y(25))))/(exp(y(7))*exp(y(7)))));
g1(19,25)=(-(exp(y(5))*params(5)*exp(y(25))/exp(y(7))));
g1(19,27)=exp(y(27));
g1(20,4)=T(23);
g1(20,8)=(-((-(exp(y(8))*exp(y(4))*(1-params(4))*exp(y(24))))/(exp(y(8))*exp(y(8)))));
g1(20,24)=T(23);
g1(20,28)=exp(y(28));
g1(21,5)=T(24);
g1(21,9)=(-((-(exp(y(9))*exp(y(5))*(1-params(5))*exp(y(25))))/(exp(y(9))*exp(y(9)))));
g1(21,25)=T(24);
g1(21,29)=exp(y(29));
g1(22,8)=(-(exp(y(8))*params(10)*exp(y(28))));
g1(22,10)=exp(y(10))*exp(y(24));
g1(22,24)=exp(y(10))*exp(y(24));
g1(22,28)=(-(exp(y(8))*params(10)*exp(y(28))));
g1(23,8)=(-(exp(y(8))*params(11)*exp(y(28))));
g1(23,12)=exp(y(12))*exp(y(25));
g1(23,25)=exp(y(12))*exp(y(25));
g1(23,28)=(-(exp(y(8))*params(11)*exp(y(28))));
g1(24,9)=(-(exp(y(9))*params(12)*exp(y(29))));
g1(24,11)=exp(y(11))*exp(y(24));
g1(24,24)=exp(y(11))*exp(y(24));
g1(24,29)=(-(exp(y(9))*params(12)*exp(y(29))));
g1(25,9)=(-(exp(y(9))*params(13)*exp(y(29))));
g1(25,13)=exp(y(13))*exp(y(25));
g1(25,25)=exp(y(13))*exp(y(25));
g1(25,29)=(-(exp(y(9))*params(13)*exp(y(29))));
g1(26,16)=(-(exp(y(16))*params(14)*exp(y(30))));
g1(26,18)=exp(y(18))*exp(y(24));
g1(26,24)=exp(y(18))*exp(y(24));
g1(26,30)=(-(exp(y(16))*params(14)*exp(y(30))));
g1(27,16)=(-(exp(y(16))*params(15)*exp(y(30))));
g1(27,20)=exp(y(20))*exp(y(25));
g1(27,25)=exp(y(20))*exp(y(25));
g1(27,30)=(-(exp(y(16))*params(15)*exp(y(30))));
g1(28,17)=(-(exp(y(17))*params(16)*exp(y(31))));
g1(28,19)=exp(y(19))*exp(y(24));
g1(28,24)=exp(y(19))*exp(y(24));
g1(28,31)=(-(exp(y(17))*params(16)*exp(y(31))));
g1(29,17)=(-(exp(y(17))*params(17)*exp(y(31))));
g1(29,21)=exp(y(21))*exp(y(25));
g1(29,25)=exp(y(21))*exp(y(25));
g1(29,31)=(-(exp(y(17))*params(17)*exp(y(31))));
g1(30,6)=T(25);
g1(30,14)=(-(params(1)*(-(exp(y(14))*exp(y(6))*params(6)*exp(y(26))))/(exp(y(14))*exp(y(14)))));
g1(30,26)=T(25);
g1(30,30)=exp(y(30))-params(1)*(1-params(8))*exp(y(30));
g1(31,7)=T(26);
g1(31,15)=(-(params(1)*(-(exp(y(15))*exp(y(7))*params(7)*exp(y(27))))/(exp(y(15))*exp(y(15)))));
g1(31,27)=T(26);
g1(31,31)=exp(y(31))-params(1)*(1-params(9))*exp(y(31));
g1(32,32)=1-params(18);
g1(33,33)=1-params(19);
g1(34,34)=1-params(20);
g1(35,35)=1-params(21);

end
