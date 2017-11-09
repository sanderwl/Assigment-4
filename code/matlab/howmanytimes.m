addpath('C:\Users\sande\OneDrive\Dokumenter\GitHub\Assigment-4\code\build-project4-Desktop_Qt_5_9_2_MSVC2017_64bit-Debug') %Change to build path
E = load('energy1.000000.txt');
E2 = load('energy2.400000.txt');

edges = unique(E);
counts = histc(E(:), edges);

edgesT24 = unique(E2);
countsT24 = histc(E2(:), edgesT24);


%%
figure(1) %T = 1.0
plot(edges,counts)
grid on
axis equal
xlabel('Energy value')
ylabel('Number of appearances')
title('How many times an energy value appears with 1 000 000 MC cycles with T = 1.0 and random initial matrix')

figure(2) %T = 2.4
plot(edgesT24,countsT24)
grid on
xlabel('Energy value')
ylabel('Number of appearances')
title('How many times an energy value appears with 1 000 000 MC cycles with T = 2.4 and random initial matrix')

%%
variance = var(E2)
