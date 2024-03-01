function[AC, NSDE, wAll] = PM(Hb1, Hd1, HbE1, HdE1, HbDesired1, method, eh, ed, belta)
%% 
if nargin<9
    belta  = 0;
end
    wAll = zeros(size(Hb1, 2), size(Hb1, 3));
    AC = zeros(1, size(Hb1, 3));
    NSDE = zeros(1, size(Hb1, 3));
    N = size(Hb1, 2);

    for i = 1:size(Hb1, 3)
        w0 = 1; % initialize virtual source filter gain(coefficient).
        pdB = HbDesired1(:, i) * w0; % desired sound pressure in bright zone.
        scaling = 1;%norm(pdB);
%         pdD = zeros(size(Hd1, 1), 1);
%         pd = [pdB; pdD];
%         Hbe = squeeze(HbE1(:, :, i))/scaling;
%         Hde = squeeze(HdE1(:, :, i))/scaling;
%         Hbd = [squeeze(Hb1(:, :, i)); squeeze(Hd1(:, :, i))]/scaling;
%         pdB = pdB/scaling;
        switch method
            case 0
                pdD = zeros(size(Hd1, 1), 1);
                pd = [pdB; pdD];
                Hbe = squeeze(HbE1(:, :, i));
                Hde = squeeze(HdE1(:, :, i));
                Hbd = [squeeze(Hb1(:, :, i)); squeeze(Hd1(:, :, i))];
                % compute Filter coefficient with Inevitable solution method
                R = (Hbd') * Hbd;
                belta = max(eig(R))*(1e-5);
                w = pinv(R+belta*eye(size(Hb1, 2))) * (Hbd') * pd;
                wAll(:, i) = w;
                PerformanceIndex = CaculateAC_NSDE;
                [AC(i), NSDE(i)] = PerformanceIndex.AC_NSDE(w, pdB, Hbe, Hde);
            case 1
                %
                Hbe = squeeze(HbE1(:, :, i))/scaling;
                Hde = squeeze(HdE1(:, :, i))/scaling;
                Hb = squeeze(Hb1(:, :, i))/scaling;
                Hd = squeeze(Hd1(:, :, i))/scaling;
                pdB = pdB/scaling;
                cvx_begin
                    variable w(N) complex
                    minimize(square_pos(norm(Hb * w - pdB)))
                    subject to
                        square_pos(norm(Hd * w)) <= ed;
                        square_pos(norm(w)) <= eh;
                cvx_end
                wAll(:, i) = w;
                PerformanceIndex = CaculateAC_NSDE;
                [AC(i), NSDE(i)] = PerformanceIndex.AC_NSDE(w, pdB, Hbe, Hde);
             otherwise
                disp('error');
                break;
        end
    end
end

