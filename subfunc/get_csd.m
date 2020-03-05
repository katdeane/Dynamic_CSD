function output = get_csd(x,dim,distan,methchoi,conduc)
% this function computes the csd of the matrix x; it is assumed that the
% values are in Volts [V]; the dist-matrix is a matrix or a scalar
% determining the intersamplingpointdistance/s - it is assumed that 
% distances are values in µm; methchoi is specifying how to compute the csd
% either: 1 - distance between points is the same (dist(1,1) is taken) or 2 
% - distance between points varies (values in dist are used accordingly)
%
% for future versions: cond is a conductivity matrix; for now set to 1

% check input
if nargin < 1 
    disp('Not enough input arguments! Program terminates!');
    return
end

if ~exist('dim')
    dim = 1;
end

if ~exist('distan') | isempty(distan)
    disp('Sampling distance not specified! Set to 0.1 mm.');
    distan = 0.1;
end

if ~exist('methchoi')
    disp('CSD is computed with equidistant points.');
    methchoi = 1;
end

% initialize
conduc = 1;


% create matrix x1, where every second column refers to the n-th (n = width-1) sampling
% point from the first column
x1        = shiftdim(x,dim-1);
x1dims    = size(x1);

if methchoi == 1
    x2 = conduc*[-diff(x1,2,1)]/(distan(1,1)^2);
elseif methchoi == 2
    if size(distan,1) ~= size(x1,1)
        disp('Size of distance matrix and size of matrix x do not match!.');
        disp('Sampling distance was set to 0.1 mm.');
        distan = [0; ones(size(x1,1)-1,1)./10];
    else
        distan = distan./1000;
    end   
    % create matrix of sum of 2 subsequent distances
    mat = distan(sort(repmat([2:size(distan)]',2,1)));mat(2:2:end) = mat(2:2:end)*-1;
    dist2 = diff(mat);dist2=dist2(2:2:end);
    
    x3 = -diff(x1,1,1)./repmat(distan(2:end),[1 x1dims(2:end)]);x4 =-diff(x3,1,1);
    x4dims = size(x4);
    x2 = ((2*conduc)./(repmat(dist2,[1 x4dims(2:end)]))).*(x4);
end

% shift dimensions back
output        = shiftdim(x2,size(x1dims,2)-dim+1);