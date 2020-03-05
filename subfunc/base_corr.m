function output = base_corr(x,basesp,dim)
% this function corrects for baseline differences; x is matrix to correct;
% basesp are the number of sampling points to use for baseline; dim is the
% dimension to work on

% check input
if nargin < 2 
    disp('Not enough input arguments! Program terminates!');
    return
end

if ~exist('dim')
    dim = 1;
end


% create matrix x1, where every second column refers to the n-th (n = width-1) sampling
% point from the first column
x1        = shiftdim(x,dim-1);
x1dims    = size(x1);

% main computation
x2 = repmat(mean(x1(1:basesp,:,:,:),1),[x1dims(1),1]);
x3 = x1-x2;

% shift dimensions back
if size(size(x3)) == size(size(x))
    output        = shiftdim(x3,size(x1dims,2)-dim+1);
else
    output        = shiftdim(x3,-(dim-1));
end

