%addpath('C:\Users\sande\OneDrive\Dokumenter\GitHub\Assigment-4\code\build-project4-Desktop_Qt_5_9_2_MSVC2017_64bit-Debug') %Change to build path

%Temperature vector
Temp = load('temp.txt');

%Files for lattice size 40*40 to 100x100
%Energy
E40 = load('energy40.txt');
E60 = load('energy60.txt');
E80 = load('energy80.txt');
E100 = load('energy100.txt');
%Magnetization
M40 = load('magnetic40.txt');
M60 = load('magnetic60.txt');
M80 = load('magnetic80.txt');
M100 = load('magnetic100.txt');
%Specific heat
Heat40 = load('heat40.txt');
Heat60 = load('heat60.txt');
Heat80 = load('heat80.txt');
Heat100 = load('heat100.txt');
%Susceptibility
Sus40 = load('sus40.txt');
Sus60 = load('sus60.txt');
Sus80 = load('sus80.txt');
Sus100 = load('sus100.txt');

%%

%Energy plots
figure(1)
hold on
plot(Temp, -E40)
plot(Temp, -E60)
plot(Temp, -E80)
plot(Temp, -E100)
legend('40x40', '60x60', '80x80', '100x100')
grid on
xlabel('Temperature')
ylabel('Energy')
title('Energy versus temperature for different lattice sizes and 10^6 MC cycles')

%Magnetization plots
figure(2)
hold on
plot(Temp, M40)
plot(Temp, M60)
plot(Temp, M80)
plot(Temp, M100)
legend('40x40', '60x60', '80x80', '100x100')
grid on
xlabel('Temperature')
ylabel('Magnetization')
title('Magnetization versus temperature for different lattice sizes and 10^6 MC cycles')

%Specific heat plots
figure(3)
hold on
plot(Temp, Heat40)
plot(Temp, Heat60)
plot(Temp, Heat80)
plot(Temp, Heat100)
legend('40x40', '60x60', '80x80', '100x100')
grid on
xlabel('Temperature')
ylabel('Specific heat')
title('Specific heat versus temperature for different lattice sizes and 10^6 MC cycles')

%Susceptibility plots
figure(4)
hold on
plot(Temp, Sus40)
plot(Temp, Sus60)
plot(Temp, Sus80)
plot(Temp, Sus100)
legend('40x40', '60x60', '80x80', '100x100')
grid on
xlabel('Temperature')
ylabel('Susceptibility')
title('Susceptibility versus temperature for different lattice sizes and 10^6 MC cycles')






