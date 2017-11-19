addpath('C:\Users\sande\OneDrive\Dokumenter\GitHub\Assigment-4\code\build-project4-Desktop_Qt_5_9_2_MSVC2017_64bit-Debug') %Change to build path

M = load('magnetization1.000000.txt'); %T = 1.0
MU = load('magnetizationUP1.000000.txt'); %T = 1.0
mcM = 100*linspace(1,length(M),length(M));
E = load('energy1.000000.txt'); %T = 1.0
EU = load('energyUP1.000000.txt'); %T = 1.0
mcE = 100*linspace(1,length(E),length(E));
M2 = load('magnetization2.400000.txt'); %T = 2.4
M2U = load('magnetizationUP2.400000.txt'); %T = 2.4
mcM2 = 100*linspace(1,length(M2),length(M2));
E2 = load('energy2.400000.txt'); %T = 2.4
E2U = load('energyUP2.400000.txt'); %T = 2.4
mcE2 = 100*linspace(1,length(E2),length(E2));

%%

%Energy plots
figure(1) %T = 1.0
hold on
plot(mcE,E)
plot(mcE2,EU)
legend('Random','Up spin','fontsize',28)
grid on
xlabel('number of mc cycles')
ylabel('Energy')
title('Energy of lattice for T = 1.0')

figure(2) %T = 2.4
hold on
plot(mcE,E2)
plot(mcE2,E2U)
legend('Random','Up spin','fontsize',28)
grid on
xlabel('number of mc cycles')
ylabel('Energy')
title('Energy of lattice for T = 2.4')

%Absolute magnetization Plotter
figure(3) %T = 1.0
hold on
plot(mcM,M)
plot(mcM,MU)
legend('Random','Up spin','fontsize',28)
ylim([399,401])
grid on
xlabel('number of mc cycles')
ylabel('Magnetization')
title('Absolute value of magnetization for T = 1.0')

figure(4) %T = 2.4
hold on
plot(mcM2,M2)
plot(mcM2,M2U)
legend('Random','Up spin','fontsize',28)
grid on
xlabel('number of mc cycles')
ylabel('Magnetization')
title('Absolute value of magnetization for T = 2.4')



