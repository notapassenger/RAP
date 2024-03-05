%% 
close all;clear;clc
addpath('functions');
addpath('toolbox');
%% Array Setup
load("array.mat");
% room
L = [4.5, 4.5, 2.2]; 
% LoundspeakerPosition
s = array.s;
% MicrophonePosition
bCtrPtsPositions = array.bCtrPtsPositions; % Bright
dCtrPtsPositions = array.dCtrPtsPositions; % Dark

% Position grapha
figure(1)
scatter(s(:, 1), s(:, 2), 10, 'filled');
hold on;
scatter(bCtrPtsPositions(:, 1), bCtrPtsPositions(:, 2), 5, 'filled');
hold on;
scatter(dCtrPtsPositions(:, 1), dCtrPtsPositions(:, 2), 5, 'filled');
hold off;
title('ArrayPosition');
grid on;
axis equal;
ATF.array = array;
%% FFT set 
c = 343;
fs = 48000;
fsResampled = fs/3;
beta = 0.3; % RT
nDFT = 3200; % n, DFT
n = 2967*3; % rir length in Time domain
LVirtualSource = 8;
%% Microphone array Rotation
% theta = 13/180;
% Rotation = [cos(theta) -sin(theta); sin(theta) cos(theta)]; 
% % bright
% bCtrPtsRot = bCtrPtsPositions(:, 1:2) - bCtrPtsPositions(19, 1:2);
% bCtrPtsNew = (Rotation * bCtrPtsRot' + bCtrPtsPositions(19, 1:2)')';
% bCtrPtsNew(:, 3) = 1.2;
% bCtrPtsPositions = bCtrPtsNew;
% % dark
% dCtrPtsRot = dCtrPtsPositions(:, 1:2) - dCtrPtsPositions(19, 1:2);
% dCtrPtsNew = (Rotation * dCtrPtsRot' + dCtrPtsPositions(19, 1:2)')';
% dCtrPtsNew(:, 3) = 1.2;
% dCtrPtsPositions = dCtrPtsNew;
%% Room and ArrayPosition(Place the array in the room)
% plot_new(a, s, bCtrPtsPositions, dCtrPtsPositions);
% strFigureSet = ["1"; "scatter"; ];
figure(2)
scatter(s(:, 1), s(:, 2), 10, 'filled');
hold on;
scatter(bCtrPtsPositions(:, 1), bCtrPtsPositions(:, 2), 5, 'filled');
hold on;
scatter(dCtrPtsPositions(:, 1), dCtrPtsPositions(:, 2), 5, 'filled');
hold on;
rectangle('Position',[-L(1) 0 L(1) L(2)]);
hold off;
title('Room and ArrayPosition');
grid on;
axis equal;
%% generate desired rir matrix, only Reserved direct path, Resampled, and FFT
directPath = zeros(1, size(bCtrPtsPositions, 1));
HB = [];
RIR = [];
for T = 1:size(bCtrPtsPositions, 1)
    distance = norm(bCtrPtsPositions(T, 1:2) - s(LVirtualSource, 1:2));
    directPath(T) = floor(distance/c * fs);
end
for T = 1:size(bCtrPtsPositions, 1)
    rir = rir_generator(c, fs, bCtrPtsPositions(T, :), s(LVirtualSource, :), L, beta, n);
    [HB(T, :), f,  RIR(T, :)] = resampled_fft(rir, nDFT, fsResampled, fs, directPath(T));                          
end
ATF.irDesired.HB = HB;
ATF.irDesired.rir = RIR;
%% generate true rir matrix, resampled(48k->16k) and FFT
HB = zeros(size(bCtrPtsPositions, 1), size(s, 1), nDFT/2 +1);
HD = HB;
HBRirResampled = zeros(size(bCtrPtsPositions, 1), size(s, 1), n/(fs/fsResampled));
HDRirResampled = HBRirResampled;
% control points
for T = 1:size(bCtrPtsPositions, 1)
    for TT = 1:size(s, 1)
        rir = rir_generator(c, fs, bCtrPtsPositions(T, :), s(TT, :), L, beta, n);
        [HB(T, TT, :), f, HBRirResampled(T, TT, :)] = resampled_fft(rir, nDFT, fsResampled, fs);
    end                           
end
for T = 1:size(dCtrPtsPositions, 1)
    for TT = 1:size(s, 1)
        rir = rir_generator(c, fs, dCtrPtsPositions(T, :), s(TT, :), L, beta, n);
        [HD(T, TT, :), f, HDRirResampled(T, TT, :)] = resampled_fft(rir, nDFT, fsResampled, fs);
    end                           
end
ATF.irTrue.HB = HB;
ATF.irTrue.HD = HD;
ATF.irTrue.HBRir = HBRirResampled;
ATF.irTrue.HDRir = HDRirResampled;
ATF.room.c = c;
ATF.room.L = L;
ATF.room.beta = beta;
ATF.irTrue.f = f;
ATF.irTrue.fs = fs;
ATF.irTrue.fsResampled = fsResampled;
ATF.irTrue.nDFT = nDFT;
%% add Posion pertubation, get Posion pertubation matrix.
DisturbanceTimesTh = 30;
HB = [];
RirResampled = [];
    % input: origin x-coordinate, origin y-coordinate, uniformly
    % distributed range, add Posion pertubation times;
    % output: Posion pertubation Coordinate;
[delta_x, delta_y] = PositionDisturbance(0, 0, 0.2, DisturbanceTimesTh); 
for i = 1:DisturbanceTimesTh
    bCtrPtsPositions = array.bCtrPtsPositions;
    bCtrPtsPositions(:, 1:2) = bCtrPtsPositions(:, 1:2) + [delta_x(i), delta_y(i)]; 
    % control points, bright zone
    for T = 1:size(bCtrPtsPositions, 1)
        for TT = 1:size(s, 1)
            rir = rir_generator(c, fs, bCtrPtsPositions(T, :), s(TT, :), L, beta, n);
            [HB(T, TT, :, i), f, RirResampled(T, TT, :, i)] = resampled_fft(rir, nDFT, fsResampled, fs);
        end                           
    end
end
% % Position grapha
% figure(1)
% scatter(s(:, 1), s(:, 2), 10, 'filled');
% hold on;
% scatter(bCtrPtsPositions(:, 1), bCtrPtsPositions(:, 2), 5, 'filled');
% hold on;
% scatter(dCtrPtsPositions(:, 1), dCtrPtsPositions(:, 2), 5, 'filled');
% hold off;
% title('ArrayPosition');
% grid on;
% axis equal;
ATF.irMeasured.HB = HB;
ATF.irMeasured.rir = RirResampled;
% save ATF.mat ATF;
%%
save ATF.mat ATF;
