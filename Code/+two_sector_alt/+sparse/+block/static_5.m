function [y, T] = static_5(y, x, params, sparse_rowval, sparse_colval, sparse_colptr, T)
  y(8)=y(13)*(y(14)/(1-params(4)))^(1-params(4));
  y(21)=y(26)*(y(27)/(1-params(10)))^(1-params(10));
end
