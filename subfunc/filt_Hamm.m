function output = filt_Hamm(x,width,dim)
% FILT_HAMM_00: this function function filters a n-dim matrix on the specified dimension
% and with the specified kernel width
%
% Syntax: output = filt_Hamm(x,width,dim);
%
% See also filt_FFT_06.

% check input
if nargin < 1 
    disp('Not enough input arguments! Program terminates!');
    return
end

if ~exist('width','var')
    width = 3;
elseif exist('width','var') && rem(width,2) == 0
    width = width+1;
    fprintf('Width is even. Width was set to: %d.\n', width);
end

if ~exist('dim','var')
    dim = 1;
end

% create matrix x3, where every second column refers to the n-th (n = width-1) sampling
% point from the first column
x1        = shiftdim(x,dim-1);
x1dims    = size(x1);

% preallocating hamm for speed
hamm = nan(size(x1,1)-(width-1),size(x1,2));
j = 1;

hammmat = repmat(hamming(width),[1,size(x1,2),size(x1,3)]);       
for i = ((width-1)/2+1):(size(x1,1)-(width-1)/2)
    hamm(j,:) = mean((x1(i-(width-1)/2:(i+(width-1)/2),:,:).*hammmat),1);
    j = j + 1;
end

% shift dimensions back
output        = shiftdim(hamm,size(x1dims,2)-dim+1);