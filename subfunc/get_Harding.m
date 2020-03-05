function [out] = get_Harding(x,dim)

% check input
if nargin < 1 
    disp('Not enough input arguments! Program terminates!');
    return
end

if ~exist('dim','var')
    dim = 1;
end

% create shifted matrix x1, where dim 1 is the working dimension
x1        = shiftdim(x,dim-1);
x1dims    = size(x1);

% main computation
x2 = sum(x1);              % computes absolute residues 
x3 = mean(x1);             % computes averaged absolute residues;
x4 = x2./sum(abs(x1));     % computes relative residues
x5 = mean(abs(x1));        % computes rectified mean (see Schroeder et al. Cereb Cortex epub ahead of print September 2006)
x6 = x1./repmat(sum(abs(x1)),[size(x1,1),1,1]); % computes relative residues for each channel individually; this enables layer specific analysis
x7 = x1./repmat(x2,[size(x1,1),1,1]); % contribution of channels to the relative residues

% shift dimensions back
out.var1 = shiftdim(x2,size(x1dims,2)-dim+1);
out.var2 = shiftdim(x3,size(x1dims,2)-dim+1);
out.var3 = shiftdim(x4,size(x1dims,2)-dim+1);       
out.var4 = shiftdim(x5,size(x1dims,2)-dim+1);
out.var5 = shiftdim(x6,size(x1dims,2)-dim+1);
out.var6 = shiftdim(x7,size(x1dims,2)-dim+1);
