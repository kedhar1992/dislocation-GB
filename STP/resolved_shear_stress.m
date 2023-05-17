% Copyright 2013 Max-Planck-Institut f�r Eisenforschung GmbH
function tau = resolved_shear_stress(eulers, n, d, sigma_a, varargin)
%% Function used to calculate the resolved shear stress
% From C.N. Reid," Deformation Geometry for Materials Scientists, Pergamon Press, Oxford, United Kingdom, 1973 (p.115-133).
% eulers: Bunge Euler angles in degrees
% d: slip direction
% n: slip plane normal
% sigma: stress in sample coordinate system
% author: c.zambaldi@mpie.de

tol = 1e-9; % tolerance for tests and for checking if normals and directions are perpendicular;
% 
% if nargin < 5
%     testFlag = false;
% end
% if nargin < 4
%     sigma = unitstress(3); % Tension along Z
% end
% if nargin < 2
%     d = [-1,1,0];
%     n = [1,1,1];
% end
% if nargin < 1
%     eulers = randBunges();
%     testFlag = true;
% end

% n = n/norm(n);
% d = d/norm(d);
% 
% if dot(n,d) > tol
%     warning_commwin('n,b not perpendicular');
% end
% 
% g = eulers2g(eulers)';
% S = schmidmatrix(d,n);
% 
% schmid = g * S * g';
% 
% % Resolved shear stress
% tau = 0;
% for ii = 1:3
%     for jj = 1:3
%         tau = tau + sigma(ii,jj) * schmid(jj,ii);
%     end
% end
% 
% 
% if testFlag
%     fprintf('Euler angles = [%.3f %.3f %.3f]\n', eulers);
%     fprintf('Resolved Shear Stress = %.4f\n', tau);
% end




%%%% Formualtions from Dr. Rajeev Kapoor

n = n/norm(n);
d = d/norm(d);
g = eulers2g(eulers);

sigma_rot = g * (sigma_a * transpose(g));
tau = abs(dot((sigma_rot * n'),d));    % n - slip plane normal;  d - burger vector 


% Verification
% tau = 0;
% for i = 1:3
%    for j = 1:3
%        tau = tau + (d(i) + sigma_a(i,j) + n(j));
%    end
% end


end 