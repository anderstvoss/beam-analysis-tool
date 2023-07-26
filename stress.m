function [Tmax, Smax, Pmax] = stress(d, v, m)

%==========================================================================
%
%   Inputs
%	    d:cm	scalar of diameter
%
%	    v:N	    scalar of shear
%
%       m:N*m	scalar of moment
%
%   Outputs
%	    v:Pa	transverse shear stress
%
%       m:Pa	bending stress
%
%       p:Pa	maximum principal stress in section
%
%   Matlab needs to be able to call pStress.m
%
%==========================================================================

%close all;

%==========================================================================
% Check input for errors
%==========================================================================

if ~isequal(class(d), 'double')
    error('Invalid diameter: diameter be a double value.');
elseif min(d) <=0
    error('Invalid diameter: diameter must be a postive nonzero value.');
end

if ~isequal(class(v), 'double')
    error('Invalid shear force: shear force must be a double value.');
end

if ~isequal(class(m), 'double')
    error('Invalid bending moment: bending moment must be a double value.');
end

%if numel(d) ~= 1 || numel(v) ~= 1 || numel(m) ~= 1
%    error('Invalid inputs: inputs must be of length 1.')
%end

%==========================================================================
% Convert to SI units
%==========================================================================

d = d * 1E-2;   % convert from cm to m

%==========================================================================
% Calculate area A and area moment of inertia I
%   - Cross section is circular
%==========================================================================

I = pi .* d .^4 ./ 64;
A = pi .* d .^2 ./ 4;

%==========================================================================
% Calculating max shear stress
%
% Relevant Formulae:
%   - tau = V*B/It
%   - B = yc * a' 
%
%   - Centroid of a segment of a circle as a function of theta
%   -- y_avg = 4 * R * sin^3(1/2 * theta) / (3 * (theta - sin(theta))
%
%   - Area of a segment of a circle as a function of theta
%   -- A = 1/2 * R^2 * (theta - sin(theta))
%
%   - Segment angle as a function of y: 0 <= y <= r
%   --theta = 2 * arccos(y/R)
%
%   - Interface area as a function of theta
%   -  t = 2 * R * sin(1/2 * theta)
%
% Assumptions:
%   - Cordinate system centered at centroid
%   - Max shear stress located at y = 0 for the circular cross section
%   -- Therefore theta = pi
%==========================================================================

Q = ( d./2 ) .* ( 4/(3*pi) ) .* ( A./2 );
B = d;

Tmax = ( v .* Q) ./ ( I .* B );

%==========================================================================
% Calculating max bending stress
%
% Relevant Formulae
%   - sigma = -M * y / I
%   - sigma_max = M * c / I
% 
% Assumptions: 
%   - c = R = d/2
%==========================================================================

Smax = ( m .* d ) ./ ( 2 .* I );

%==========================================================================
% Calculating max principal stress
%   - use pStress(d, v, m);
%==========================================================================

Pmax = pStress(d, v, m);