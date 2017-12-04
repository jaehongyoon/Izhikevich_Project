clc; clear; % clear cmd window and cache
addpath('./func') % add path of funtions
if ~exist('./results')
    mkdir('./results')
end
    
for I = 0.3:0.1:1.5
    bifurcation_information = IK_run_ab_coordinates(I);
    sorted_data = IK_draw_plot(bifurcation_information);
    
    save(sprintf('./results/bif_data_I_%f.mat', I), 'bifurcation_information','sorted_data');
    pause();
end