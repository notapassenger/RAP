clear all; clc
%%
load ATF.mat;
addpath('functions');
addpath('algorithm');
%% 
% get array size.
Mb = size(ATF.array.bCtrPtsPositions, 1); % the number of microphones or called control points in bright zone.
Md = size(ATF.array.dCtrPtsPositions, 1);
N = size(ATF.array.s, 1); % the number of Loudspeakers array.  

% get Hb and Hd.
Fs = ATF.irTrueResampled.fsResampled;% sample frequency
f = ATF.irTrueResampled.f; % true Analog frequency
fSpace = 1; % removed frequency bins, User variable
f_max = floor(length(f)/fSpace);% The total number of frequency points taken
HB = ATF.irTrueResampled.HB;
HD = ATF.irTrueResampled.HD;
HBDesired = ATF.HBDesiredResampled;
% HBDesired = squeeze(ATF.irTrueResampled.HB(:, 8, :));
HBMeasured = ATF.irMeasured.HB;
NoisyNumber = size(HBMeasured, 4);
% HBE = ATF.irTrue.HBE;
% HDE = ATF.irTrue.HDE;
% LVirtualSource = 8; % 8th

% The ir value of the frequency point was extracted
Hb = HB(:, :, 1:fSpace:end); % control
% HbE = HBE(:, :, 1:fSpace:end); % evaluation
Hd = HD(:, :, 1:fSpace:end);
% HdE = HDE(:, :, 1:fSpace:end);
HbMeasured = HBMeasured(:, :, 1:fSpace:end, :); % noisy
% HbDesired = squeeze(HB(:, LVirtualSource, 1:fSpace:end)); % desired pressure ir
HbDesired = HBDesired(:, 1:fSpace:end); % desired pressure ir
clear HB HBE HD HDE HBMeasured HB;


%% AC
% AC_ACCNoisyPre = zeros(size(Hb, 3), NoisyNumber);
[AC_ACC] = ACC(Hb, Hd, Hb, Hd, 0);
% [AC_ACCN] = ACC(Hb, Hd, HbE, HdE, 1);
% for i = 1:NoisyNumber
%     [AC_ACCNoisyPre(:, i)] = ACC(squeeze(HbMeasured(:, :, :, i)), Hd, HbE, HdE, 0);
% end
% AC_ACCNoisy = mean(AC_ACCNoisyPre, 2);
%%
figure(1)
plot(f(1:fSpace:end), AC_ACC);
% hold on;
% plot(f(1:fSpace:end), AC_ACCNoisy);
% hold on;
% plot(f(1:fSpace:end), AC_RACC);
% hold on;
% plot(f(1:fSpace:end), AC_ACCN);
% hold off;
% xlabel('f/Hz');
% ylabel('dB');
% legend('True', 'Measured', 'Opt', 'N');
grid on;
xlim([100 4000]);
set(gca,'XScale','log');
title('AC');

%% PM
% eh = db2mag(13); % the power of virtual source
eh = db2mag(13)^2;
% ed = db2mag(-30);
ed = db2mag(-21)^2;
NoisyNumber = 1;
% startnew = 11;
% endnew = 20;
Index = 1:1:801;
% AC_PMNoisyPre = zeros(size(HbMeasured, 3), NoisyNumber);
Lindex = length(Index);
AC_PMNoisyPre = zeros(NoisyNumber, Lindex);
NSDE_PMNoisyPre = zeros(NoisyNumber, Lindex);
[AC_PM, NSDE_PM, w_PM] = PM(Hb(:, :, Index), Hd(:, :, Index), ...
    Hb(:, :, Index), Hd(:, :, Index), HbDesired(:, Index), 0, eh, ed);
%%
for i = 1:NoisyNumber
    [AC_PMNoisyPre(i, :), NSDE_PMNoisyPre(i, :), ~] = PM(squeeze(HbMeasured(:, :, Index, i+19)), ...
        Hd(:, :, Index), Hb(:, :, Index), Hd(:, :, Index), HbDesired(:, Index), 0, eh, ed);
end
AC_PMNoisy = mean(AC_PMNoisyPre, 1);
NSDE_PMNoisy = mean(NSDE_PMNoisyPre, 1);
%%
% Robust PM
% eh = db2mag(13); % the power of virtual source
% ed = db2mag(-30);
[AC_RPMPre, NSDE_RPMPre, ~] = RobustPM(Hb(:, :, Index), HbMeasured(:, :, Index, 20:(19+NoisyNumber)), ...
    Hd(:, :, Index), Hb(:, :, Index), Hd(:, :, Index), HbDesired(:, Index), eh, ed);
AC_RPM = mean(AC_RPMPre, 1);
NSDE_RPM = mean(NSDE_RPMPre, 1);
%% F norm Test
% meanHbMeasured = mean(HbMeasured(:, :, :, 20:25), 4);
% deltaHb = meanHbMeasured - Hb;
% normDeltaHb = 0;
% for i = 1:size(meanHbMeasured, 3)
%     normDeltaHb = normDeltaHb + norm(squeeze(deltaHb(:, :, i)), 'fro');
% end
%%
figure(8)
plot(f(Index), AC_PM(1:Lindex));
% hold on;
% plot(f(Index), AC_PMNoisy(1:Lindex));
% hold on;
% plot(f(Index), AC_RPM(1:Lindex));
hold off;
xlabel('f/Hz');
ylabel('dB');
% legend('True', 'Measured', 'Opt');
% legend('True', 'Measured');
grid on;
xlim([100 4000]);
set(gca,'XScale','log');
title('AC');

figure(9)
plot(f(Index), NSDE_PM(1:Lindex));
% hold on;
% plot(f(Index), NSDE_PMNoisy(1:Lindex));
% hold on;
% plot(f(Index), NSDE_RPM(1:Lindex));
hold off;
xlabel('f/Hz');
ylabel('dB');
% legend('True', 'Measured');
% legend('True', 'Measured', 'Opt');
grid on;
xlim([100 4000]);
set(gca,'XScale','log');
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
%%
figure(2)
% plot(f(1:fSpace:end), mean(AC_RPM, 1));
plot(f(1:50), AC_PM(1:50));
hold on
plot(f(1:50), AC_PMNoisyC(1:50));
hold on
plot(f(1:50), AC_RPM);
hold off;
xlabel('f/Hz');
ylabel('dB');
legend('true', 'measured', 'Opt');
grid on;
xlim([100 4000]);
set(gca,'XScale','log');
title('AC');

figure(3)
% plot(f(1:fSpace:end), mean(NSDE_RPM, 2));
plot(f(1:50), NSDE_PM(1:50));
hold on
plot(f(1:50), NSDE_PMNoisyC(1:50));
hold on
plot(f(1:50), NSDE_RPM);
hold off;
xlabel('f/Hz');
ylabel('dB');
legend('true', 'measured', 'Opt');
grid on;
xlim([100 4000]);
set(gca,'XScale','log');
title('NSDE');
%% Robust ACC-PM


