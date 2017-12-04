function sorted_info = IK_draw_plot(bif_struct)
% sort the structure object accordingly to type of bifurcation
sorted_info = sort_bif(bif_struct);
color_code = ['b','r','g','k','m','c','y'];
str = {};
figure(100)
clf;
hold on
for i = 1:numel(sorted_info)
    plot(sorted_info(i).a, sorted_info(i).b, 'Color', color_code(i))
    str{end+1} = sorted_info(i).type;
end
hold off
xlabel('a');ylabel('b');legend(str, 'Location', 'Southwest')

end

function bif = sort_bif(bif_struct)
% sort the bifurcation information object accordingly to the fxpt type
% and bifurcation
bif = [];
for i = 1:numel(bif_struct) % check all objects in info
    info = bif_struct(i);
    type = bifurcation_type(info);
    
    [TF, idx] = in_structure_bif(bif, type);
    if TF % if the type of bifurcation exist as field of bif
        bif(idx).a(end+1) = info.a_b(1);
        bif(idx).b(end+1) = info.a_b(2);
    else
        bif(end+1).type = type;
        bif(end).a = info.a_b(1);
        bif(end).b = info.a_b(2);
    end
    
end
end

function type = bifurcation_type(info)
% return the bifurcation type
fxpt_type_current = info.fxpt_type_current; % the current fxpt type
fxpt_type_previous = info.fxpt_type_previous; % the previous fxpt type

if isempty(fxpt_type_current) || isempty(fxpt_type_previous)
    % if either is empty (no fxpt), its a saddle point bifurcation
    type = 'saddle node bifurcation';
else
    stability_current = cell(1, numel(fxpt_type_current));
    stability_previous = cell(1, numel(fxpt_type_previous));
    
    % do stability test
    for i = 1:numel(fxpt_type_previous)
        stability_current{i} = stability(fxpt_type_current{i});
    end
    
    for i = 1:numel(fxpt_type_previous)
        stability_previous{i} = stability(fxpt_type_previous{i});
    end
    
    % supercritical pitchfork bif
    if (numel(fxpt_type_previous) == 1 && numel(fxpt_type_current) == 3 && ...
            isempty(setdiff(stability_current, {'stable','unstable','stable'}))...
            && isempty(setdiff(stability_previous, {'stable'}))) ||... 
            (numel(fxpt_type_current) == 1 && numel(fxpt_type_previous) == 3 && ...
            isempty(setdiff(fxpt_type_previous, {'stable','unstable','stable'}))...
            && isempty(setdiff(fxpt_type_current, {'stable'})))
        type = 'supercritical pitchfork';
    % subcritical pitchfork bif
    elseif (numel(fxpt_type_previous) == 1 && numel(fxpt_type_current) == 3 && ...
            isempty(setdiff(stability_current, {'unstable','stable','unstable'}))...
            && isempty(setdiff(stability_previous, {'unstable'}))) ||... 
            (numel(fxpt_type_current) == 1 && numel(fxpt_type_previous) == 3 && ...
            isempty(setdiff(fxpt_type_previous, {'unstable','stable','unstable'}))...
            && isempty(setdiff(fxpt_type_current, {'unstable'})))
        type = 'subcritical pitchfork';
    % supertranscritical bif
    elseif isempty(setdiff(stability_current, flip(stability_previous))) &&...
            ~strcmp(stability_current{1}, stability_current{2}) && ...
            numel(stability_current) == 2
        type = 'supertranscritical bif';
    % subtranscritical bif
    elseif isempty(setdiff(stability_current, stability_previous)) &&...
            ~strcmp(stability_current{1}, stability_current{2}) && ...
            numel(stability_current) == 2
        type = 'subtranscritical bif';
    else
        % homoclinic bif
        % same number of fxpt but the fxpt stability changes
        if numel(stability_current) == numel(stability_previous) && ...
                ~isempty(setdiff(stability_current, stability_previous)) 
            type = 'homoclinic bif';
        else % unknowns
            type = 'unknown';
        end
        
    end
    
end
end

function s = stability(str)
if strcmp(str, 'stable node') || strcmp(str, 'stable spirals')
    s = 'stable';
else
    s = 'unstable';
end
end

function [TF, idx] = in_structure_bif(bif, type)
% structure comparison function for fixed points

% check if current data in fxpt array
for i = 1:length(bif)
    if( strcmp(bif(i).type, type) )        
        % if true, return values
        idx = i;
        TF = true;
        return        
    end
end

TF = false;
idx = NaN;

end