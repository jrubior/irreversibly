function [y, T] = dynamic_5(y, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T)
  y(43)=y(48)*(y(49)/(1-params(4)))^(1-params(4));
  y(56)=y(61)*(y(62)/(1-params(10)))^(1-params(10));
end
