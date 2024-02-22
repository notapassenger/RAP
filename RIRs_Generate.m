%% 
close all;clear;clc

%% Array Setup
% room
L = [4.9, 7.22, 5.29]; 
% loudspeakersPosition
l = 1:30;
R = 1.5;
xL = -1 * R * sin((2 * l - 1) * pi / 30);
yL = R * cos((2 * l - 1) * pi / 30);

% MicrophonePosition
r1 = 0.12;
r2 = 0.10;
m = 1:60; % 30 control points, 30 Evaluation points
rm = 0.12-0.02*floor((m-1)/30);
% A zone, Bright
xmA = -1*rm.*sin(mod(m-1, 30)*2*pi/30)-0.4;
ymA = rm.*cos(mod(m-1, 30)*2*pi/30);
% B zone, Dark
xmB = -1*rm.*sin(mod(m-1, 30)*2*pi/30)+0.4;
ymB = rm.*cos(mod(m-1, 30)*2*pi/30);
% Control points
xmAControl = xmA(1:2:end);
ymAControl = ymA(1:2:end);
xmBControl = xmB(1:2:end);
ymBControl = ymB(1:2:end);
% Evaluation points
xmAEvaluation = xmA(2:2:end);
ymAEvaluation = ymA(2:2:end);
xmBEvaluation = xmB(2:2:end);
ymBEvaluation = ymB(2:2:end);

% Position grapha
figure(1)
scatter(xL, yL, 10, 'filled');
hold on;
scatter(xmAControl, ymAControl, 5, 'filled');
hold on;
scatter(xmBControl, ymBControl, 5, 'filled');
hold on;
scatter(xmAEvaluation, ymAEvaluation, 5, 'filled');
hold on;
scatter(xmBEvaluation, ymBEvaluation, 5, 'filled');
hold off;
title('ArrayPosition');
grid on;
axis equal;

%% Room and ArrayPosition(Place the array in the room)
arrayCenter = [2.45, 3.8, 1.45];% Array Center
s = [xL', yL', 1.45*ones(size(xL, 2), 1)]+arrayCenter;
bCtrPtsPositions = [xmAControl', ymAControl', 1.45*ones(size(xmAControl, 2), 1)]+arrayCenter;
dCtrPtsPositions = [xmBControl', ymBControl', 1.45*ones(size(xmBControl, 2), 1)]+arrayCenter;
bEvaPtsPositions = [xmAEvaluation', ymAEvaluation', 1.45*ones(size(xmAEvaluation, 2), 1)]+arrayCenter;
dEvaPtsPositions = [xmBEvaluation', ymBEvaluation', 1.45*ones(size(xmBEvaluation, 2), 1)]+arrayCenter;
figure(2)
scatter(s(:, 1), s(:, 2), 10, 'filled');
hold on;
scatter(bCtrPtsPositions(:, 1), bCtrPtsPositions(:, 2), 5, 'filled');
hold on;
scatter(dCtrPtsPositions(:, 1), dCtrPtsPositions(:, 2), 5, 'filled');
hold on;
scatter(bEvaPtsPositions(:, 1), bEvaPtsPositions(:, 2), 5, 'filled');
hold on;
scatter(dEvaPtsPositions(:, 1), dEvaPtsPositions(:, 2), 5, 'filled');
hold on;
rectangle('Position',[0 0 L(1) L(2)]);
hold off;
title('Room and ArrayPosition');
grid on;
axis equal;
ATF.array.bCtrPtsPositions = bCtrPtsPositions;
ATF.array.dCtrPtsPositions = dCtrPtsPositions;
ATF.array.bEvaPtsPositions = bEvaPtsPositions;
ATF.array.dEvaPtsPositions = dEvaPtsPositions;
ATF.array.s = s;
                 
%% generate true rir matrix,  and convert to the frequency domain
c = 343;
fs = 48000;
beta = 0.3; % RT
n = 4096; % n, DFT
HB = [];
HD = [];
% control points
for T = 1:size(bCtrPtsPositions, 1)
    for TT = 1:size(s, 1)
        rir = rir_generator(c, fs, bCtrPtsPositions(T, :), s(TT, :), L, beta, n);
        Y = fft(rir);
        P2 = (Y/n);
        P1 = P2(1:floor(n/2) + 1);
        P1(2:end-1) = 2*P1(2:end-1);
        f = fs*(0:floor(n/2))/n; 
        HB(T, TT, :) = P1';
    end                           
end
for T = 1:size(dCtrPtsPositions, 1)
    for TT = 1:size(s, 1)
        rir = rir_generator(c, fs, dCtrPtsPositions(T, :), s(TT, :), L, beta, n);
        Y = fft(rir);
        P2 = (Y/n);
        P1 = P2(1:floor(n/2) + 1);
        P1(2:end-1) = 2*P1(2:end-1);
        f = fs*(0:floor(n/2))/n; 
        HD(T, TT, :) = P1';
    end                           
end
% evaluation points
for T = 1:size(bEvaPtsPositions, 1)
    for TT = 1:size(s, 1)
        rir = rir_generator(c, fs, bEvaPtsPositions(T, :), s(TT, :), L, beta, n);
        Y = fft(rir);
        P2 = (Y/n);
        P1 = P2(1:floor(n/2) + 1);
        P1(2:end-1) = 2*P1(2:end-1);
        f = fs*(0:floor(n/2))/n; 
        HBE(T, TT, :) = P1';
    end                           
end
for T = 1:size(dEvaPtsPositions, 1)
    for TT = 1:size(s, 1)
        rir = rir_generator(c, fs, dEvaPtsPositions(T, :), s(TT, :), L, beta, n);
        Y = fft(rir);
        P2 = (Y/n);
        P1 = P2(1:floor(n/2) + 1);
        P1(2:end-1) = 2*P1(2:end-1);
        f = fs*(0:floor(n/2))/n; 
        HDE(T, TT, :) = P1';
    end                           
end
ATF.irTrue.HB = HB;
ATF.irTrue.HD = HD;
ATF.irTrue.HBE = HBE;
ATF.irTrue.HDE = HDE;
ATF.room.c = c;
ATF.room.L = L;
ATF.room.beta = beta;
ATF.FFT.fs = fs;
ATF.FFT.n = n;
ATF.FFT.f = f;


%% generate rir matrix adding perturbation with microphone positions, and convert to the frequency domain 
%  mismatch in bright zone, 11 groups matrix Noisy.

addNoisyNumber = 10;
perturbationGrid = -0.2 + (0.2 + 0.2) .* rand(addNoisyNumber,2);
HBPertur = [];
for i = 1:addNoisyNumber
    bPerturbation = bCtrPtsPositions + [perturbationGrid(i, 1), perturbationGrid(i, 2), 0];
    for T = 1:size(bCtrPtsPositions, 1)
        for TT = 1:size(s, 1)
            rir = rir_generator(c, fs, bPerturbation(T, :), s(TT, :), L, beta, n);           
            Y = fft(rir);
            P2 = (Y/n);
            P1 = P2(1:floor(n/2) + 1);
            P1(2:end-1) = 2*P1(2:end-1);
            f = fs*(0:floor(n/2))/n; 
            HBPertur(T, TT, :, i) = P1';
        end                           
    end
end
ATF.irMeasured.HB = HBPertur;
ATF.irMeasured.addNoisyNumber = addNoisyNumber;
ATF.irMeasured.perturbationGrid = perturbationGrid;
save ATF.mat ATF;


