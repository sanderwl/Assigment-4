addpath('C:\Users\sande\OneDrive\Dokumenter\GitHub\Assigment-4\code\build-project4-Desktop_Qt_5_9_2_MSVC2017_64bit-Debug') %Change to build path
M = load('magnetizationT1.txt'); %T = 1.0
E = load('energyT1.txt'); %T = 1.0
M2 = load('magnetizationT24.txt'); %T = 2.4
E2 = load('energyT24.txt'); %T = 2.4

%%
mc = linspace(1,length(E),length(E));
mc2 = linspace(1,length(E2),length(E2));

mc3 = linspace(1,length(M),length(M));
mc4 = linspace(1,length(M2),length(M2));


%%

% figure(1)
% histogram(E, mc)
% xlabel('mc cycles')
% ylabel('Energy')

%Energy plots
figure(2) %T = 1.0
plot(mc,E)
ylim([-850,-750])
xlabel('mc cycles / 50')
ylabel('Energy')
title('Energy vs Monte Carlo cycles for T = 1.0')

figure(3) %T = 2.4
plot(mc2,E2)
xlabel('mc cycles / 50')
ylabel('Energy')
title('Energy vs Monte Carlo cycles for T = 2.4')

%Magnetization plots
figure(4) %T = 1.0
plot(mc3,-M)
ylim([350,450])
xlabel('mc cycles / 50')
ylabel('Magnetization')
title('Magnetization vs Monte Carlo cycles for T = 1.0')

figure(5)
plot(mc4,-M2)
xlabel('mc cycles / 50')
ylabel('Magnetization')
title('Magnetization vs Monte Carlo cycles for T = 2.4')




