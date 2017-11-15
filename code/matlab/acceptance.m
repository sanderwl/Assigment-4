addpath('C:\Users\sande\OneDrive\Dokumenter\GitHub\Assigment-4\code\build-project4-Desktop_Qt_5_9_2_MSVC2017_64bit-Debug') %Change to build path

Acc = load('TempAcceptance.txt');
A = load('acceptanceMC1.000000.txt');
AU = load('acceptanceMCUP1.000000.txt');
mcA = 100*linspace(1,length(A),length(A));
A2 = load('acceptanceMC2.400000.txt');
A2U = load('acceptanceMCUP2.400000.txt');
mcA2 = 100*linspace(1,length(A2),length(A2));

%Acceptance plot
figure(1)
plot(Acc(:,1),Acc(:,2))
grid on
xlabel('Temperature')
ylabel('Acceptance')
title('Acceptance for mcs = 1000 and random initial matrix')

figure(2) %T = 1.0
hold on
plot(mcA,A)
plot(mcA,AU)
legend('Random','Up spin')
grid on
xlabel('number of mc cycles')
ylabel('Acceptance')
title('Acceptance for T = 1.0 and random initial matrix')

figure(3)
hold on
plot(mcA2,A2)
plot(mcA2,A2U)
legend('Random','Up spin')
grid on
xlabel('Number of mc cycles')
ylabel('Acceptance')
title('Acceptance for T = 2.4 and random initial matrix')
