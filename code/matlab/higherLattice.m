%addpath('C:\Users\sande\OneDrive\Dokumenter\GitHub\Assigment-4\code\build-project4-Desktop_Qt_5_9_2_MSVC2017_64bit-Debug')
%The files needed are already in this folder

%Temperature vector
Temp = load('temp.txt');

%Files for lattice size 40*40 to 100x100
%Energy
E40 = load('energy40.txt')/(40*40);
E60 = load('energy60.txt')/(60*60);
E80 = load('energy80.txt')/(80*80);
E100 = load('energy100.txt')/(100*100);
%Magnetization
M40 = load('magnetic40.txt')/(40*40);   
M60 = load('magnetic60.txt')/(60*60);
M80 = load('magnetic80.txt')/(80*80);
M100 = load('magnetic100.txt')/(100*100);
%Specific heat
Heat40 = load('heat40.txt')/(40*40);
Heat60 = load('heat60.txt')/(60*60);
Heat80 = load('heat80.txt')/(80*80);
Heat100 = load('heat100.txt')/(100*100);
%Susceptibility
Sus40 = load('sus40.txt')/(40*40);
Sus60 = load('sus60.txt')/(60*60);
Sus80 = load('sus80.txt')/(80*80);
Sus100 = load('sus100.txt')/(100*100);


%%

%Energy plots
figure(1)
hold on
plot(Temp, E40)
plot(Temp, E60)
plot(Temp, E80)
plot(Temp, E100)
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

%%

%Calculating the maximum value of specific heat and susceptibility
[HM40, HI40] = max(Heat40);
[HM60, HI60] = max(Heat60);
[HM80, HI80] = max(Heat80);
[HM100, HI100] = max(Heat100);

[SM40, SI40] = max(Sus40);
[SM60, SI60] = max(Sus60);
[SM80, SI80] = max(Sus80);
[SM100, SI100] = max(Sus100);

TempArray = [Temp(HI40),Temp(HI60),Temp(HI80),Temp(HI100)];
L = [40,60,80,100];

figure(5)
plot((1./L),TempArray)
grid on
xlabel('1/L')
ylabel('Critical temperature')
title('Critical temperature against 1/L')





