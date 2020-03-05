function SSr = SS(varargin)

vec = varargin{1};
m = varargin{2};
SSr = sum((vec - m) .^ 2);
