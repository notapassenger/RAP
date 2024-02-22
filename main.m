clear all; clc
%%
load ATF.mat;
addpath('toolbox');
addpath('functions');
addpath('algorithm');
%% 
% get array size.
Mb = size(ATF.array.bCtrPtsPositions, 1); % the number of microphones or called control points in bright zone.
Md = size(ATF.array.dCtrPtsPositions, 1);
N = size(ATF.array.s, 1); % the number of Loudspeakers array.  

% get Hb and Hd.
Fs = ATF.FFT.fs;% sample frequency
f = ATF.FFT.f; % true Analog frequency
fSpace = 20; % removed frequency bins, User variable
f_max = floor(length(f)/fSpace);% The total number of frequency points taken
HB = ATF.irTrue.HB;
HD = ATF.irTrue.HD;
HBMeasured = ATF.irMeasured.HB;
NoisyNumber = ATF.irMeasured.addNoisyNumber-1;
HBE = ATF.irTrue.HBE;
HDE = ATF.irTrue.HDE;
LVirtualSource = 8; % 8th

%% The ir value of the frequency point was extracted
Hb = HB(:, :, 1:fSpace:end); % control
HbE = HBE(:, :, 1:fSpace:end); % evaluation
Hd = HD(:, :, 1:fSpace:end);
HdE = HDE(:, :, 1:fSpace:end);
HbMeasured = HBMeasured(:, :, 1:fSpace:end, :); % noisy
HbDesired = squeeze(HB(:, LVirtualSource, 1:fSpace:end)); % desired pressure ir
clear HB HBE HD HDE HBMeasured HB;


%% AC
AC_ACCNoisyPre = zeros(size(Hb, 3), NoisyNumber);
[AC_ACC] = ACC(Hb, Hd, HbE, HdE, 0);
[AC_ACCN] = ACC(Hb, Hd, HbE, HdE, 1);
for i = 1:NoisyNumber
    [AC_ACCNoisyPre(:, i)] = ACC(squeeze(HbMeasured(:, :, :, i)), Hd, HbE, HdE, 0);
end
AC_ACCNoisy = mean(AC_ACCNoisyPre, 2);
%%
figure(1)
plot(f(1:fSpace:end), AC_ACC);
hold on;
plot(f(1:fSpace:end), AC_ACCNoisy);
hold on;
plot(f(1:fSpace:end), AC_RACC);
hold on;
plot(f(1:fSpace:end), AC_ACCN);
hold off;
xlabel('f/Hz');
ylabel('dB');
legend('True', 'Measured', 'Opt', 'N');
grid on;
title('AC');

%% PM
eh = 1; % the power of virtual source
AC_PMNoisyPre = zeros(size(Hb, 3), NoisyNumber);
NSDE_PMNoisyPre = zeros(size(Hb, 3), NoisyNumber);
[AC_PM, NSDE_PM, wAllPM] = PM(Hb, Hd, HbE, HdE, HbDesired, eh);
for i = 1:NoisyNumber
    [AC_PMNoisyPre(:, i), NSDE_PMNoisyPre(:, i)] = PM(squeeze(HbMeasured(:, :, :, i)), Hd, HbE, HdE, HbDesired, eh);
end
AC_PMNoisy = mean(AC_PMNoisyPre, 2);
NSDE_PMNoisy = mean(NSDE_PMNoisyPre, 2);
%
figure(2)
plot(f(1:fSpace:end), AC_PM);
hold on;
plot(f(1:fSpace:end), AC_PMNoisy);
hold off;
xlabel('f/Hz');
ylabel('dB');
legend('True', 'Measured');
grid on;
title('AC');
figure(2)
plot(f(1:fSpace:end), NSDE_PM);
hold on;
plot(f(1:fSpace:end), NSDE_PMNoisy);
hold off;
xlabel('f/Hz');
ylabel('dB');
legend('True', 'Measured');
grid on;
title('NSDE');
%% ACC-PM
eh = 1; % the power of virtual source
karpa = 0.9; % Weight 
AC_APNoisyPre = zeros(size(Hb, 3), NoisyNumber);
NSDE_APNoisyPre = zeros(size(Hb, 3), NoisyNumber);
[AC_AP, NSDE_AP] = ACC_PM(Hb, Hd, HbE, HdE, HbDesired, eh, karpa);
for i = 1:NoisyNumber
    [AC_APNoisyPre(:, i), NSDE_APNoisyPre(:, i)] = ACC_PM(squeeze(HbMeasured(:, :, :, i)), Hd, HbE, HdE, HbDesired, eh, karpa);
end
AC_APNoisy = mean(AC_APNoisyPre, 2);
NSDE_APNoisy = mean(NSDE_APNoisyPre, 2);
%
figure(3)
plot(f(1:fSpace:end), AC_AP);
hold on;
plot(f(1:fSpace:end), AC_APNoisy);
hold off;
xlabel('f/Hz');
ylabel('dB');
legend('True', 'Measured');
grid on;
title('AC');
figure(4)
plot(f(1:fSpace:end), NSDE_AP);
hold on;
plot(f(1:fSpace:end), NSDE_APNoisy);
hold off;
xlabel('f/Hz');
ylabel('dB');
legend('True', 'Measured');
grid on;
title('NSDE');

%% Robust ACC
AC_RACCNoisyPre = zeros(size(Hb, 3), NoisyNumber);
[AC_RACC] = RobustACC(Hb, Hd, HbE, HdE);
for i = 1:NoisyNumber
    [AC_RACCNoisyPre(:, i)] = RobustACC(squeeze(HbMeasured(:, :, :, i)), Hd, HbE, HdE);
end
AC_RACCNoisy = mean(AC_RACCNoisyPre, 2);
%
figure(5)
plot(f(1:fSpace:end), AC_RACC);
hold on;
plot(f(1:fSpace:end), AC_RACCNoisy);
hold off;
xlabel('f/Hz');
ylabel('dB');
legend('True', 'Measured');
grid on;
title('AC');
%% Robust PM

%% Robust ACC-PM


