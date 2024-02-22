function[AC, NSDE] = ACC_PM(Hb1, Hd1, HbE1, HdE1, HbDesired1, eh, karpa)
    for i = 1:size(Hb1, 3)
    %     scale = 1e-4;
        w0 = 1; % initialize virtual source filter gain(coefficient).
        pdB = HbDesired1(:, i) * w0; % desired sound pressure in bright zone.
        Hbe = squeeze(HbE1(:, :, i));
        Hde = squeeze(HdE1(:, :, i));
        Hb2 = squeeze(Hb1(:, :, i));
        Hd2 = squeeze(Hd1(:, :, i));
        % compute Filter coefficient with Inevitable solution method
        R = karpa*(Hd2')*Hd2 + (1-karpa)*(Hb2')*Hb2;
        belta = max(eig(R))*(1e-5);
        w = pinv(R + belta*eye(size(Hb2, 2))) * (1-karpa) * (Hb2') * pdB;  
    %     wAll(:, i) = w;
        PerformanceIndex = CaculateAC_NSDE;
        [AC(i), NSDE(i)] = PerformanceIndex.AC_NSDE(w, pdB, Hbe, Hde, i);
    end
end
