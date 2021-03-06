addpath('C:\Users\sande\OneDrive\Dokumenter\GitHub\Assigment-4\code\build-project4-Desktop_Qt_5_9_2_MSVC2017_64bit-Debug') %Change to build path

E = load('energy1.000000.txt');
EU = load('energyUP1.000000.txt');
E2 = load('energy2.400000.txt');
E2U = load('energyUP2.400000.txt');

edges = unique(E);
edgesU = unique(EU);
counts = histc(E(:), edges);
countsU = histc(EU(:), edgesU);

edgesT24 = unique(E2);
edgesT24U = unique(E2U);
countsT24 = histc(E2(:), edgesT24);
countsT24U = histc(E2U(:), edgesT24U);

%%

figure(1) %T = 1.0
hold on
plot(edges,counts)
plot(edgesU,countsU)
grid on
xlabel('Energy value')
ylabel('Number of appearances')
title('How many times an energy value appears with 1 000 000 MC cycles with T = 1.0 and random initial matrix')

figure(2) %T = 2.4
hold on
plot(edgesT24,countsT24)
plot(edgesT24U,countsT24U)
grid on
xlabel('Energy value')
ylabel('Number of appearances')
legend('Distribution random', 'Distribution up spin')
title('How many times an energy value appears with 1 000 000 MC cycles with T = 2.4 and random initial matrix')

variance1 = var(E)
variance2 = var(E2)
