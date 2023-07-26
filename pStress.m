function [Pmax] = pStress(d, v, m, ny)

%==========================================================================
%
%   Inputs
%	    d:m	    vector of diameters
%
%	    v:N	    vector of shears
%
%	    m:N*m	vector of moments
%
%       (optional)
%	    ny	    number of y sections to evaluate principal stress at
%
%   Outputs
%	    p:Pa	maximum principal stress in each section
%
%==========================================================================

%==========================================================================
% Configuration
%==========================================================================

NY = 101;

%==========================================================================
% Check input for errors
%==========================================================================

if ~isequal(class(d), 'double')
    error('Invalid diameter: diameter vector must be of double values.');
elseif min(d) <=0
    error('Invalid diameter: diameter vector must contain only postive nonzero values.');
end

if ~isequal(class(v), 'double')
    error('Invalid shear force: shear force vector must be of double values.');
end

if ~isequal(class(m), 'double')
    error('Invalid bending moment: bending moment vector must be a double value.');
end

if length(unique([numel(d), numel(v), numel(m)])) > 1
    error('Invalid input vectors: input vectors must be of identical length.')
end
len = numel(d);

if ~exist('ny','var')
    ny = NY;
end

%==========================================================================
% Parse and clean inputs
%==========================================================================

%reshape to ensure row vector
d = reshape(d, [1, len]);
v = reshape(v, [1, len]);
m = reshape(m, [1, len]);

%==========================================================================
% Generate matrix of y diameter y values
%   -d/2 <= y <= d/2
%==========================================================================

y = zeros(ny, len);
for i = 1:len
    y(:,i) = linspace(-d(i)/2, d(i)/2, ny);
end

%==========================================================================
% Calculate area moment of inertia I
%   - Cross section is circular
%==========================================================================

I = pi .* d .^4 ./ 64;  % vector of cross sections

%==========================================================================
% Calculate Matrix of Normal and Shear Stresses
%   - sigma_xx = bending stress = -My/I
%   - sigma_xy = shear stresss = VQ/IB
%==========================================================================

% normal stress
sigma_xx = -1 .* m .* y ./ I;
%sigma_yy = 0;

% shear stress
% Q = y_avg * A
% angle of segmaent theta = 2 .* acos(Y ./ r)
% area of segment A = 1/2 * r^2 .* (theta - sin(theta))
% centroid of segment y_avg = 4 * r * (sin(1/2 * theta) .^3) ./ (3 .* (theta - sin(theta)))
% interface of segment B = 2 * r .* sin(1/2 .* theta)
% y_avg * A / B = 1/3 * r^2 * sin(arccos(y/r))^2

sigma_xy = ( v./I ) .* ( 1/3 ) .* ( d.^2 ./ 4 ) .* sin(acos( 2.*y./d )).^2;

%==========================================================================
% Calculate principal stresses
%   - p1 = (s_x+s_y)/2 +- sqrt( (s_x-s_y)^2/4 + s_xy^2 )
%==========================================================================

p1 = sigma_xx./2 + sqrt( (sigma_xx./2).^2 + sigma_xy.^2 );

%==========================================================================
% Sort for max principal stress
%==========================================================================

Pmax = max(p1, [], 1);