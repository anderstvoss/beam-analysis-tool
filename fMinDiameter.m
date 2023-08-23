%==========================================================================
% Setup
%
% Notes:
%   - This script may take a few minutes or longer to run
%   - Higher values for nx, ny, or nd exponentially increase runtime
%   - dmin, dmax, and nd need to be chosen carefully in order to get accurate
%     results 
%   - Matlab needs to be able to call minDia.m, pStress.m and vmInt.m
%==========================================================================
    
X = 15;         % cm, location of support
L = 30;         % cm, length of beam
W = -120;       % N/cm, distributed load
PMAX = 290;     % MPa, max principal stress
nx = 100;       % number of sections to analyze
ny = 101;       % number of y sections to evalate principal stress at

dmin = 0.01;    % cm, minimum diameter to evaluate
dmax = 2;       % cm, maximum diameter to evaluate
nd = 1E3;       % number of diameters to evaluate

% vector of diameters to evaluate
pD = linspace(dmin, dmax, nd);      

%==========================================================================
% Convert to SI units
%==========================================================================

x = X * 1E-2;           % convert from cm to m
l = L * 1E-2;           % convert from cm to m
w = W * 1E2;            % convert from N/cm to N/m
Pmax = PMAX * 1E6;      % convert from MPa to Pa
pD = pD .* 1E-2;        % convert from cm to m

%==========================================================================
% Calculate minimum diameters
%==========================================================================

md = minDia(x, l, w, Pmax, nx, ny, pD);     % m, vector of minimum diameters
md = md * 1E2;      % convert from m to cm
minimalDiameter = md;

%==========================================================================
% Plot minimum diameters
%==========================================================================

xp = linspace(0, L, nx);    % x coordinate of cross section, cm
z = zeros(1,nx);
padding = 0.1;              % padding for plot

fig = figure;
hold on;

plot(xp, md);
plot(xp, z, 'k');
plot(X, 0, "k*");
title(sprintf('Diameter vs X'));
subtitle(sprintf('Support at x = %d cm', X));
xlabel('Position (cm)');
ylabel('Minimum Diameter (cm)');
xlim([0,L]);
ylim([-padding, max(md)+padding]);
hold off