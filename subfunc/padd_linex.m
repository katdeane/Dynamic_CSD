function output = padd_linex(x,siz,dim)
%this function handles 3D-matrices
%
% version change: - 2010-01-17 changed direction of gradient at end of matrix
% (paddmat2) added the replicated matrix instead of subtracting it

% check input
if nargin < 1 
    disp('Not enough input arguments! Program terminates!');
    return
end

if ~exist('siz')
    siz = 3;
    disp(sprintf('Size was set to %d more sampling points on each end of the matrix.', siz));  
end

if ~exist('dim')
    dim = 1;
end

% shift dimension to the according to the working dimension
x1        = shiftdim(x,dim-1);
x1dims    = size(x1);

% padd new points
paddmat1 = repmat(x1(1,:,:),[siz,1,1])-(repmat(flipud([1:siz]'),[1 size(x1,2),size(x1,3)]).*(repmat(diff(x1([1:2],:,:)),[siz,1,1])));
paddmat2 = repmat(x1(end,:,:),[siz,1,1])+(repmat([1:siz]',[1 size(x1,2),size(x1,3)]).*(repmat(diff(x1([end-1:end],:,:)),[siz,1,1])));
paddmat = [paddmat1; x1; paddmat2];

% shift dimensions back
output        = shiftdim(paddmat,size(x1dims,2)-dim+1);