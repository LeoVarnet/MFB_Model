function [ output_args ] = display_cfg( cfg, varargin )
%DISPLAY_CFG fprintf function for configuration parameters
%   display_cfg( cfg ) or display_cfg( cfg, 'name1','name2'... )
% Displays on the sceen the values of some configuration parameters, or all
% parameters if no name is specified 
% 
% Léo Varnet 2016

if isempty(varargin)
    varargin = fieldnames(cfg);
end
for ifield = 1:length(varargin)
    fieldvalue = getfield(cfg, varargin{ifield});
    if isnumeric(fieldvalue)
        stringtoplot = num2str(fieldvalue(:)');
    elseif isstr(fieldvalue)
        stringtoplot = fieldvalue;
    else
        stringtoplot = '...';
    end
    fprintf([' - ' varargin{ifield} ': ' stringtoplot '\n']);
end

end

