addpath('C:\Users\sande\OneDrive\Dokumenter\GitHub\Assigment-4\code\build-project4-Desktop_Qt_5_9_2_MSVC2017_64bit-Debug') %Change to build path

%Loading energy arrays
E1=load('energy1.000000.txt');
E125=load('energy1.250000.txt');
E15=load('energy1.500000.txt');
E175=load('energy1.750000.txt');
E2=load('energy2.000000.txt');
E225=load('energy2.250000.txt');
E25=load('energy2.500000.txt');
E275=load('energy2.750000.txt');
E3=load('energy3.000000.txt');

%%

%Making array with energy values with no repeats
edges1 = unique(E1);
edges125 = unique(E125);
edges15 = unique(E15);
edges175 = unique(E175);
edges2 = unique(E2);
edges225 = unique(E225);
edges25 = unique(E25);
edges275 = unique(E275);
edges3 = unique(E3);

%Calculating counts per energy value
counts1 = histc(E1(:), edges1);
counts125 = histc(E125(:), edges125);
counts15 = histc(E15(:), edges15);
counts175 = histc(E175(:), edges175);
counts2 = histc(E2(:), edges2);
counts225 = histc(E225(:), edges225);
counts25 = histc(E25(:), edges25);
counts275 = histc(E275(:), edges275);
counts3 = histc(E3(:), edges3);

%Plotting the probability
figure(9)
hold all
plot (edges1,counts1,'linewidth',2);
plot (edges125,counts125,'linewidth',2);
plot (edges15,counts15,'linewidth',2);
plot (edges175,counts175,'linewidth',2);
plot (edges2,counts2,'linewidth',2);
plot (edges225,counts225,'linewidth',2);
plot (edges25,counts25,'linewidth',2);
plot (edges275,counts275,'linewidth',2);
plot (edges3,counts3,'linewidth',2);
title('Random initial matrix for varying temperature','fontsize',28)
l=legend('1','1.25','1.5','1.75','2.0','2,25','2,5','2.75','3.0')
l.FontSize=28;
ylabel('number of occurences')
xlabel('Energy')
grid on
xlim([-800,-100])