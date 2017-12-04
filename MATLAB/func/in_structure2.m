function [TF, idx] = in_structure2(fxpt, type, bf_type, a, I)
% Author: Jaehong Yoon
% Date: 09/29/17
% structure comparison function for fixed points

% check if current data in fxpt array
for i = 1:length(fxpt)
    if( (fxpt(i).type == type) && strcmp(fxpt(i).bf_type, bf_type) && ...
            (fxpt(i).a == a) && (fxpt(i).I == I))        
        % if true, return values
        idx = i;
        TF = true;
        return        
    end
end

TF = false;
idx = NaN;

end