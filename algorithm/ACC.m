function[AC] = ACC(Hb1, Hd1, HbE1, HdE1, DL)
% DL diagonal loading, 0 or 1
    AC = zeros(1, size(Hb1, 3));
    for i = 1:size(Hb1, 3)
    %     scale = 1e-4;
        Hbe = squeeze(HbE1(:, :, i));
        Hde = squeeze(HdE1(:, :, i));
        % compute Filter coefficient with Inevitable solution method
        Rb = (squeeze(Hb1(:, :, i))') * squeeze(Hb1(:, :, i));
        Rd = (squeeze(Hd1(:, :, i))') * squeeze(Hd1(:, :, i));
        if DL == 1
%             [EigvalueB, vectorB] = MaxEigenvector(Rb);
            EigvalueB = max(eig(Rd));
            gammaB = 0.5 * EigvalueB;
            gammaD = 0.05 * norm(Rd)/1e+2;
        else
            gammaB = 0;
            gammaD = 0;
        end
        I = eye(size(Rd));
        belta = max(eig(Rd))*(1e-5);
        R = pinv(Rd + belta*eye(size(Hd1, 2)) + gammaD*I) * (Rb - gammaB*I);
        [~, w] = MaxEigenvector(R);
    %     wAll(:, i) = w;
    %     AC(i) = 10*log10(MaxMu);
        PerformanceIndex = CaculateAC_NSDE;
        [AC(i)] = PerformanceIndex.AC(w, Hbe, Hde, i);
    end
end


