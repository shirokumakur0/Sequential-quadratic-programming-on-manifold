clear;
% This is the test case for makeCanonicalBasis.m, which makes the canonical
% basis according to the dimension of an inout manifold. Here we confirm
% the case that the input is spheres, stiefels, and obliques.

% NOTICE: we don't assume that the dimension is not zero, necessarly
% In particular, the sphere and oblique tests don't assume that the input 
% dimension equals to zero because this is a just prototyping.
% If you want to modify this, it might help you read the code of
% spherefactory at first because the spherefactory code might not asuume
% the case that the input equals zero, already. (and in that case, we can't
% do nothing unless modify the original manopt code.)

% Test 1 sphere case
for n = 1:3:30
    manifold = spherefactory(n);
    problem.M = manifold;
    tanspcdim = numel(problem.M.zerovec());
    canonicalbasis = makeCanonicalBasis(problem);
    assert(isa(canonicalbasis, 'cell'),'Sphere: The data type is not cell.')
    for idx = 1:tanspcdim
        idxbase = zeros(tanspcdim,1);
        idxbase(idx) = 1;
        assert(isequal(canonicalbasis{idx},idxbase),... 
        'Sphere: canonicalbase does not match the expected one: n: %d, idx: %d.',...
        n, idx);
    end
end

% Test 2 Stiefel manifold case
for n = 0:3:10
    for p = 0:2:n-1
        manifold = stiefelfactory(n,p);
        problem.M = manifold;
        tanspcdim = numel(problem.M.zerovec());
        canonicalbasis = makeCanonicalBasis(problem);
        assert(isa(canonicalbasis, 'cell'),'Stiefel: The data type is not cell.');
        if n == 0 || p == 0
            test = cell(n,0);
            assert(isequal(canonicalbasis, test),...
                'Stiefel: canonicalbase does not match the expected one: n: %d, idx: %d.',...
                n, 0);
        else
            for idx = 1:tanspcdim
                idxbase = problem.M.zerovec();
                idxbase(idx) = 1;
                assert(isequal(canonicalbasis{idx},idxbase),... 
                'Stiefel: canonicalbase does not match the expected one: n: %d, idx: %d.',...
                n, idx);
            end
        end
    end    
end

% Test 3 Oblique manifold case
for n = 1:3:10
    for m = 1:2:10
        manifold = obliquefactory(n,m);
        problem.M = manifold;
        tanspcdim = numel(problem.M.zerovec());
        canonicalbasis = makeCanonicalBasis(problem);
        assert(isa(canonicalbasis, 'cell'),'Oblique: The data type is not cell.');
        for idx = 1:tanspcdim
            idxbase = problem.M.zerovec();
            idxbase(idx) = 1;
            assert(isequal(canonicalbasis{idx},idxbase),... 
            'Oblique: canonicalbase does not match the expected one: n: %d, idx: %d.',...
            n, idx);
        end
    end
end

fprintf('All tests have been accepted! [test_makeCanonicalBasis]\n')