clear;
% This is a test script for setLocalDefaults.m, which set some basic
% hyperparameters for SQP algorithm. Since the function to be tested is 
% very simple and its return values might be changed when we tune the
% hyperparameters, we only test here whether the return components are
% generated desirably, or not.

% Test case 1:
problem= exampleProblemSphere();
localdefaults = setLocalDefaults(problem);
assert(isfield(localdefaults, "maxouteriter"), 'maxouteriter does not exists in the localdefaults.')
assert(isfield(localdefaults, "maxtime"), 'maxtime does not exists in the localdefaults.')
assert(isfield(localdefaults, "minstepsize"), 'minstepsize does not exists in the localdefaults.')
assert(isfield(localdefaults, "tolgradnorm"), 'tolgradnorm does not exists in the localdefaults.')
assert(isfield(localdefaults, "storedepth"), 'storedepth does not exists in the localdefaults.')
assert(isfield(localdefaults, "tau"), 'tau does not exists in the localdefaults.')
assert(isfield(localdefaults, "rho"), 'rho does not exists in the localdefaults.')
assert(isfield(localdefaults, "beta"), 'beta does not exists in the localdefaults.')
assert(isfield(localdefaults, "gamma"), 'gamma does not exists in the localdefaults.')
assert(isfield(localdefaults, "mus"), 'mus does not exists in the localdefaults.')
assert(isfield(localdefaults, "lambdas"), 'lambdas does not exists in the localdefaults.')
assert(isfield(localdefaults, "ls_max_steps"), 'ls_max_steps does not exists in the localdefaults.')

fprintf('All tests have been accepted! [test_setLocalDefaults]')