function [minDiam] = minDia(x, l, w, P, nx, ny, d1, d2, nd)

%==========================================================================
%
%   Inputs
%	    x:m	    scalar of support location
%	
%       l:m	    scalar of beam length
%	    
%       w:N/m	scalar of distributed load
%	
%       P:Pa	scalar of maximum principal stress
%	
%       (optional)
%       nx	    number of x sections to evaluate minimum diameter at
%
%       ny	    number of y sections to evaluate principal stress at
%	
%       d1:m	vector of diameters to evaluate stress of
%           OR
%       d1:m    minimum diameter to evalutate
%
%       d2:m    maximum diameter to evaluate
%
%       nd      number of diameters to evaluate
%       
%   Outputs
%	    minD:m	vector of minimum diameters
%
%   Matlab needs to be able to call pStress.m and vmInt.m
%
%==========================================================================

%==========================================================================
% Configuration
%==========================================================================

NX = 100;
NY = 101;

D1 = 0.001;
D2 = 1;
ND = 100000;

%==========================================================================
% Check input for errors
%==========================================================================

if ~isequal(class(x), 'double')
    error('Invalid support location: scalar must be a double value.');
elseif min(x) < 0
    error('Invalid support location: must not be negative.');
elseif max(x) > l
    error('Invalid support location: must not longer than the length.');
end

if ~isequal(class(l), 'double')
    error('Invalid length: length scalar must be a double value.');
elseif l <= 0
    error('Invalid length: length scalar must be a postive nonzero value.');
end

if ~isequal(class(w), 'double')
    error('Invalid dist load: distributed load scalar must be a double value.');
end

if ~isequal(class(P), 'double')
    error('Invalid stress: max stress scalar must be a double value.');
elseif P <= 0
    error('Invalid stress: max stress scalar must be a postive nonzero value.');
end

if ~exist('nx','var')
    nx = NX;
end

if ~exist('ny','var')
    ny = NY;
end

if exist('d1','var') && numel(d1) > 1
    d = d1;
    nd = numel(d1);
elseif exist('d1','var') && exist('d2','var') && exist('nd','var')
    d = linspace(d1, d2, nd);
elseif  exist('d1','var') && exist('d2','var')
    d = linspace(d1, d2, ND);
    nd = ND;
else
    d = linspace(D1, D2, ND);
    nd = ND;
end

%==========================================================================
% Parse and clean inputs
%==========================================================================

% reshape d to ensure col vector
d = reshape(d, [nd, 1]); 

%==========================================================================
% Calculate internal forces
%==========================================================================

[v, m] = vmInt(x, l, w, nx);

%==========================================================================
% Calculate principal stresses
%==========================================================================

p1 = zeros(nd, nx);

for i = 1:nd
    D = ones(1, nx) .* d(i);
    p1(i,:) = pStress(D, v, m, ny);
end

%==========================================================================
% Find diameter closest to max stress
%==========================================================================

[~, ind] = min(abs(p1 - P), [], 1);
minDiam = d(ind);