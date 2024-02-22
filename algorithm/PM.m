function[AC, NSDE, wAll] = PM(Hb1, Hd1, HbE1, HdE1, HbDesired1, eh)
%% 
    wAll = zeros(size(Hb1, 2), size(Hb1, 3));
    AC = zeros(1, size(Hb1, 3));
    NSDE = zeros(1, size(Hb1, 3));
    for i = 1:size(Hb1, 3)
    %     scale = 1e-4;
        w0 = 1; % initialize virtual source filter gain(coefficient).
        pdB = HbDesired1(:, i) * w0; % desired sound pressure in bright zone.
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
        [AC(i), NSDE(i)] = PerformanceIndex.AC_NSDE(w, pdB, Hbe, Hde, i);

    %     % compute Filter coefficient with Cvx method
    %     cvx_begin
    %         variable w(N) complex
    %         minimize(norm(Hbd*w - pd))
    %         subject to
    %               norm(w) <= eh;
    %     cvx_end
    % %     wcvxall(:, i) = w;
    %     [AC(i), NSDE(i)] = CaculateAC_NSDE(w, pdB, HbE, HdE, i);
    end
end

