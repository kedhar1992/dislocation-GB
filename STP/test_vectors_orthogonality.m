% Copyright 2013 Max-Planck-Institut f�r Eisenforschung GmbH
function res = test_vectors_orthogonality(v1, v2, tol)
%% Check orthogonality between 2 vectors (v1 and v2)

if nargin < 3
    tol = 1e-9; % tolerance for tests and for checking if normals and directions are perpendicular;
end

if nargin < 1
    for ii = 1 : 1e4
        v1 = random_direction() * rand;
        v2 = orthogonal_vector(v1);
        %         v2 = random_direction();
        %         assert(~test_vectors_orthogonality(v1,v2));
        res = test_vectors_orthogonality(v1,v2);
    end
    return
end

assert(numel(v1) == 3)
assert(numel(v2) == 3)

res = (dot(v1,v2) < tol);

end