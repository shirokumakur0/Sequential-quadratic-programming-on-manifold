function params = set_default_para(params, name, default_value, type, a, b)
% To set parameters in the struct "params".
% If the given field does not exist, then set the field to be the default value, otherwise, check whether the given value in the required range.
% 
% INPUT:
% params :a struct that need to be checked
% name : the name of the given field
% dafault_value : the default value of the field
% type : it is either "float" or "int". It is to indicate the type of the field
% [a, b] : indicate the required range of the given value.
%
% By Wen Huang
    if(~isfield(params, name))
        params = setfield(params, name, default_value);
    end
    value = getfield(params, name);

    if(strcmp(type, 'int') && round(value) ~= value)
        msg = sprintf('Invalid arguments: params.%s must be a integer', name);
        error(msg);
    end
    
    if(value < a || value > b)
        msg = sprintf('Invalid arguments: params.%s must be in [%d, %d]', name, a, b);
        error(msg);
    end
end
