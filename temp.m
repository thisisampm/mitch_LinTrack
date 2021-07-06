function temp(data, varargin)

varargin

for i=1:length(varargin)
    if ischar(varargin{i})
        varargin{i}
        if strcmp(varargin{i}, 'dance')
            disp('success')
        end
    end
end
            