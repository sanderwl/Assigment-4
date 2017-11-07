E = load('energy.txt');
mc = length(E);

%%

figure(1)
histogram(E, mc)
xlabel('mc cycles')
ylabel('Energy')

figure(2)
plot(mc,E)
xlabel('mc cycles')
ylabel('Energy')