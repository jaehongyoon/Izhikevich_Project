%%  directly find bifurcation point
I_list = 0.:0.02:0.4;
figure; hold on; box on;
for i=1:length(I_list)
    I = I_list(i);
    a = [0:0.0001:0.45]; % a samples 
    v = (a-5)/0.08;     % v samples so that tau = 0
    b = (I+140+5*v+0.04*v.^2)./v;
    if I > 0.2
        plot(a,b,'k-')
    else
        plot(a,b,'r-')
    end
%     pause;
end
ylabel('b')
xlabel('a')
%%  overlay genetic algorithm results
A = csvread('gene_A_record.csv');
B = csvread('gene_B_record.csv');
score = csvread('gene_score_record.csv');
scatter(A,B,10,bsxfun(@times,[1,0,0],(1-(score)/40.)),'filled')
% scatter(A,B,score,'.','markerfacecolor',(1-(score)/40.))
