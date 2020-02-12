function check_function_handle(fns, name)
% Check whether the required function handle exists or not.
% if it does not exist, an error message is given.
% INPUT:
% fns : a struct that contains function handles.
% name : the name of the required function handle.
% 
% By Wen Huang
    if(~isfield(fns, name) || ~isa(fcnchk(getfield(fns, name)), 'function_handle'))
        msg = sprintf('Invalid arguments: missing fns.%s or fns.%s is not a function handle', name, name);
        error(msg);
    end
end
