function sorted = bifurcation_animation(fxpt, filename)
% draw the animation of how bifurcation changes accordingly to the change
% of external current and change in a. The bifurcation diagram of b vs. v
% of Izhikevich model is drawn

sorted = sort_fxpt(fxpt); % sort data

h = figure(100);
xlabel('b');ylabel('v')
for i = 1:numel(sorted) % number of I
    for j = 1:numel(sorted(i).data) % number of a
        leg = {};
        hold on
        for k = 1:numel(sorted(i).data(j).graph_data)
            c = color_code(sorted(i).data(j).graph_data(k).type);
            plot(sorted(i).data(j).graph_data(k).b, ...
                sorted(i).data(j).graph_data(k).v, 'Color', c)
            leg{end+1} = sorted(i).data(j).graph_data(k).type;
        end
        hold off
        title(sprintf('phase diagram for I = %12.5f, a = %12.5f',sorted(i).I, sorted(i).data(j).a))
        legend(leg)
        
        % capture the plot as image
        frame = getframe(h);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256); 
        
        % Write to the GIF File 
        if i*j == 1 
            imwrite(imind, cm, filename, 'gif', 'Loopcount', inf); 
        else 
            imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append'); 
        end 
    end
end


end

function sorted = sort_fxpt(fxpt)
% sort the given fxpt object accordingly to hierarchy suited for graphics.
% The data are sorted accordingly to I in first layer, where it is sorted
% accordingly to the value of b inside it.

% sorted object (struct): first layer
% |--------.I:  external currenet
% |--------.data:   2nd layer
%           |-------.: value of a
%           |-------.graph_data: third layer
%                   |-------.b: array of b
%                   |-------.v: array of fxpt value of v accordingly to the
%                               b in same order
%                   |-------.type: fxpt type

sorted = []; % initialize sorted obj
h = waitbar(0, 'sorting data');

for i = 1:numel(fxpt) % construct the sorted structure
    waitbar(i/numel(fxpt));
    I = fxpt(i).I;
       
    idx_I = find_indx_I(sorted, I); % index of I in 1st layer of sorted obj
    if isempty(idx_I) % not found
        sorted(end+1).I = I; sorted(end).data = [];
        idx_I = numel(sorted);
    end
    
        
    idx_a = find_indx_a(sorted(idx_I), fxpt(i).a); % index of a in 2nd layer of sorted obj
    if isempty(idx_a) % not found
        sorted(idx_I).data(end+1).a = fxpt(i).a;
        sorted(idx_I).data(end).graph_data = [];
        idx_a = numel(sorted(idx_I).data);
    end
    
    sorted(idx_I).data(idx_a).graph_data(end+1).b = fxpt(i).b;
    sorted(idx_I).data(idx_a).graph_data(end).v = fxpt(i).v_star;
    sorted(idx_I).data(idx_a).graph_data(end).type = fxpt(i).bf_type;
end

close(h);

end


function idx = find_indx_I(sorted, I)
idx = [];
for i = 1:numel(sorted)
    if sorted(i).I == I
        idx = i;
        return
    end
end
end

function idx = find_indx_a(sorted, a)
idx = [];
for i = 1:numel(sorted.data)
    if sorted.data(i).a == a
        idx = i;
        return
    end
end
end

function c = color_code(str)
    if strcmp(str, 'saddle point')
        c = 'r';
    elseif strcmp(str, 'unstable node')
        c = 'g';
    elseif strcmp(str, 'saddle point')
        c = 'b';
    elseif strcmp(str, 'stable node')
        c = 'y';
    elseif strcmp(str, 'non-isolated fixed node')
        c = 'm';
    elseif strcmp(str, 'center')
        c = 'k';
    elseif strcmp(str, 'unstable spirals')
        c = 'm';
    elseif strcmp(str, 'stable spirals')
        c = 'c';
    else % homoclinic
        c = 'k';

    end
end