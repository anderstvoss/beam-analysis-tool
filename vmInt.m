function [v, m] = vmInt(x, l, w, nx, caseNoSupport)

%==========================================================================
%
%   Inputs
%	    x:m	    vector of support locations
%
%       l:m	    scalar of beam length
%	    
%       w:N/m	scalar of distributed load (upward is postive)
%	    
%       (optional)
%       nx:	    number of x sections to evaluate internal stress and moment at
%
%       caseNoSupport	    treat cases where the suppot is at 0 or l as no
%                           support
%   
%   Outputs
%	    v:N	    matrix of internal shear forces
%	
%       m:N*m	matrix of internal bending moments
%
%==========================================================================

%==========================================================================
% Configuration
%==========================================================================

NX = 100;
CNS = false;

%==========================================================================
% Check input for errors
%==========================================================================

if ~isequal(class(x), 'double')
    error('Invalid support location: vector must be of double values.');
elseif min(x) < 0
    error('Invalid support location: must not contain negative values.');
elseif max(x) > l
    error('Invalid support location: must not contain values longer than the given length.');
end

if ~isequal(class(l), 'double')
    error('Invalid length: length scalar must be a double value.');
elseif l <= 0
    error('Invalid length: length scalar must be a postive nonzero value.');
end

if ~isequal(class(w), 'double')
    error('Invalid dist load: distributed load scalar must be a double value.');
end

if ~exist('n','var')
    nx = NX;
end
if ~exist('caseNoSupport','var')
    caseNoSupport = CNS;
end

%==========================================================================
% Parse and clean inputs
%==========================================================================

% reshape x to ensure col vector
cases = numel(x);
x = reshape(x, [cases, 1]);

% generate x locations of sections
xc = linspace(0, l, nx);     

%==========================================================================
% Determine reaction force at middle support
% 
% Deflection assumptions:
% Assume superposition
% Postive forces act in the upwards y direction
%
% Simply supported beam, with distrbuted load w across entire length L
% v = -wx/24EI * (x^3 - 2L*x^2 + L^3)
%   0 <= x <= L
%
% Simply supported beam, with concentrated load P at distance a of length L
% v = -Pbx/6EIL * (L^2 - b^2 - x^2)
%   0 <= x <= a
%
% v(a) = 0
% v(a) = -wa/24EI * (a^3 - 2L*a^2 + L^3) + -Pba/6EIL * (L^2 - b^2 - a^2)
% 0 = -w/24 * (a^3 - 2L*a^2 + L^3) - Pb/6L * (L^2 - b^2 - a^2)
% Pb/6L * (L^2 - b^2 - a^2) = -w/24 * (a^3 - 2L*a^2 + L^3)
% P = -6wL/24b * (a^3 - 2L*a^2 + L^3) / (L^2 - b^2 - a^2)
%==========================================================================

a = x;
b = l - a;

P = (-6*w*l) ./ (24.*b) .* (a.^3 - 2*l .* a.^2 + l^3) ./ (l^2 - b.^2 - a.^2);

%==========================================================================
% Calculate reaction forces
%   Net force on the entire beam
%       0 = ( w * l ) + P + R_A + R_B
%   Net moment about A
%       0 = x*P + l*R_B + ( 1/2 * w *l^2 )
%==========================================================================

R_B = ( (-1/2 * w * l^2) - x.*P ) / l;
R_A = (-1*w*l) - P - R_B;

%==========================================================================
% Calculate internal shear forces
%   Net internal force on a particular section 0 <= xc < x
%       0 = R_A - V + w*xc
%       V = R_A + w*xc
%
%   Net internal force on a particular section x <= xc <= l
%       0 = R_A + P + w*xc - V
%       V = R_A + P + w*xc
%==========================================================================

v = (R_A + w .* xc ).*(xc < x) + (R_A + P + w .* xc ).*(xc >= x);

%==========================================================================
% Calculate internal bending moment
%   Net internal moment on a particular section 0 <= xc < x
%       0 = M + ( 1/2 * w * xc^2 ) - xc*V
%       M = ( -1/2 * w * xc^2 ) + xc*( R_A + w*xc )
%
%   Net internal moment on a particular section x <= xc <= l
%       0 = M + x*P + ( 1/2 * w * xc^2 ) - xc*V
%       M = ( -1/2 * w * xc^2 ) + xc*( R_A + P + w*xc) - x*P
%==========================================================================


m = (( -1/2 .* w .* xc.^2 ) + xc.*( R_A + w .* xc)).*(xc < x) + ...
    (( -1/2 .* w .* xc.^2 ) + xc.*( R_A + P + w .* xc ) - x.*P).*(xc >= x);

%==========================================================================
% Override for cases where x = 0 or x = l
%   treat as simply supported
%
%   Net force on the beam
%       0 = R_A + R_B + w*l
%   Net moment on the beam
%       0 = ( 1/2 * w * l^2 ) + l* R_B
%
%   R_A = R_B = ( -1/2 * w * l )
%
%   Net internal force on the beam
%       0 = R_A - V + w*xc
%       V = R_A + w*xc
%
%   Net internal moment on the beam
%       0 = M + ( 1/2 * w * xc^2 ) - xc*V
%       M = ( -1/2 * w * xc^2 ) + xc*( R_A + w*xc )
%==========================================================================

%find index where x = 0 or l
if caseNoSupport
    inds = [find(x == 0), find(x == l)];
    if ~isempty(inds)
        R = ( -1/2 * w * l );
        V = R + xc.*w;
        M = ( -1/2 * w * xc.^2 ) + xc.*( R + w .* xc );
        for i = inds
            v(i,:) = V;
            m(i,:) = M;
        end 
    end
end
