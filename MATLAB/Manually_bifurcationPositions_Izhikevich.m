%%  directly find bifurcation point
I_list = 0.2:0.02:0.4;
figure; hold on; box on;
for i=1:length(I_list)
    I = I_list(i);
    a = [0:0.0001:0.45]; % a samples 
    v = (a-5)/0.08;     % v samples so that tau = 0
    b = (I+140+5*v+0.04*v.^2)./v;
    plot(a,b,'k-')
    pause;
end
ylabel('b')
xlabel('a')
