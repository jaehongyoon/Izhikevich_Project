%%  directly find bifurcation point
I_list = 0.:0.02:0.4;
figure; hold on; box on;
%   PART1: hopf-bifurcation
for i=1:length(I_list)
    I = I_list(i);
    a = [0.001:0.001:0.45]; % a samples 
    v = (a-5)/0.08;     % v samples so that tau = 0
    b = (I+140+5*v+0.04*v.^2)./v;
    delta = -(.08*v+5).*a+a.*b;
    plot(a(delta>0),b(delta>0),'k-')
%     pause;
end
%   PART2: saddle-node bifurcation
I_list = 0.:0.02:0.4;
for i=1:length(I_list)
    I = I_list(i);
    a = [0.001:0.001:0.45]; % a samples 
    v0 = zeros(size(a));
    fz = @(v,aa,II) 0.04*v.^2+v.*(0.08*v+5-aa).^2./(4*aa)-140-II;
    for j = 1:length(a)
        v0(j) = fzero(@(v)fz(v,a(j),I),-59);
    end
    b = 0.08*v0+5+(0.08*v0+5-a).^2./(4*a);
    plot(a(delta<0),b(delta<0),'r-')
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
