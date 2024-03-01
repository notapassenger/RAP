function[AC, NSDE, wAll] = RobustPM(HbTrue1, HbMeasured1, Hd1, HbE1, HdE1, HbDesired1, eh, ed)
%% 
    wAll = zeros(size(HbMeasured1, 2), size(HbMeasured1, 3));
    AC = zeros(size(HbMeasured1, 4), size(HbMeasured1, 3));
    NSDE = zeros(size(HbMeasured1, 4), size(HbMeasured1, 3));
    N = size(HbMeasured1, 2);
    M = size(HbMeasured1, 1);
    iterations = 100;
%     HbstopDifferenceValue = 1e-8;
%     wstopDifferenceValue = 1e-4;
%     meanHbMeasured = squeeze(mean(HbMeasured1, 4));
    for i = 1:size(HbMeasured1, 3)
        w0 = 1; % initialize virtual source filter gain(coefficient).
        pdB = HbDesired1(:, i) * w0; % desired sound pressure in bright zone.
        scaling = 1;%norm(pdB);
        pdB = pdB/scaling;
        Hd = squeeze(Hd1(:, :, i))/scaling;
        HbTrue = squeeze(HbTrue1(:, :, i))/scaling;
        Hbe = squeeze(HbE1(:, :, i))/scaling;
        Hde = squeeze(HdE1(:, :, i))/scaling;
%% get C, ellipsoid        
%         C_pre = reshape(squeeze(HbMeasured1(:, :, i, :)), M*N, size(HbMeasured1, 4));
%         C = cov(C_pre');
%         ppC = chol(C);
%         pe = norm(ppC * reshape(HbTrue - squeeze(meanHbMeasured(:, :, i)), M*N, 1));        
%% get C, sphere
        C = eye(M*N);
        choles = chol(C);
%         pe = norm(choles \ reshape(HbTrue - meanHbMeasured(:, :, i)/scaling, M*N, 1));
        %ppC = inv(chol(C));
%% opt 
        for k = 1:size(HbMeasured1, 4)
            HbMeasured = squeeze(HbMeasured1(:, :, i, k))/scaling;
            pe = norm(choles \ reshape(HbTrue - HbMeasured, M*N, 1));
            cvx_begin
                variable wOpt(N) complex 
                minimize(square_pos(norm(HbMeasured * wOpt - pdB)))
                subject to
                    square_pos(norm(Hd * wOpt)) <= ed;
                    square_pos(norm(wOpt)) <= eh;
            cvx_end
%           wTS(u, :) = wOpt;
            %wTS(1, :) = wOpt;
            wCur = wOpt;
        % calculate G1 using w1.
            cvx_begin
                variable HbOpt(M, N) complex
                minimize(norm(HbOpt * wOpt - pdB))
                subject to
                    norm(choles \ reshape(HbOpt - HbMeasured, M * N, 1)) <= pe;
            cvx_end
            %HbTS(:, :, 1) = HbOpt;
            HbCur = HbOpt;
            for u = 2:iterations    
                % Hbd_opti(:, 1) = HTS(u - 1, :);
                cvx_begin
                    variable wOpt(N) complex
                    minimize(square_pos(norm(HbOpt * wOpt - pdB)))
                    subject to
                        square_pos(norm(Hd * wOpt)) <= ed;
                        square_pos(norm(wOpt)) <= eh;
                cvx_end
                %wTS(u, :) = wOpt;
%                     wOpt(:, 1) = wTS(u - 1, :);
                wPre = wCur;
                wCur = wOpt;
                wOpt = wPre;
                cvx_begin
                    variable HbOpt(M, N) complex
                    minimize(norm( HbOpt * wOpt - pdB))
                    subject to
                        norm(choles \ reshape(HbOpt - HbMeasured, M * N, 1)) <= pe;
                cvx_end
                %HbTS(:, :, u) = HbOpt;
                HbPre = HbCur;
                HbCur = HbOpt;
                HbOpt = HbPre;
                % Hbd_opti_all(:, :, u, T) = Hbd_opti;
%                     PerformanceIndex = CaculateAC_NSDE;
%                     [AC(u, i), NSDE(u, i)] = PerformanceIndex.AC_NSDE(wOpt, pdB, Hbe, Hde, i);              
%                 if(mod(u-2, 3) == 0 && norm(wCur-wPre)<=wstopDifferenceValue)    
%                     break;
%                 end
%                 if(mod(u, 3) == 0 && norm(HbCur-HbPre, 'fro')<=HbstopDifferenceValue)
%                     break;
%                 end
            end
            PerformanceIndex = CaculateAC_NSDE;
            [AC(k, i), NSDE(k, i)] = PerformanceIndex.AC_NSDE(wOpt, pdB, Hbe, Hde);
            % wAll(:, i) = wOpt;
        end
    end
end

