function[AC] = RobustACC(Hb1, Hd1, HbE1, HdE1)
    N = size(Hb1, 2);
    AC = zeros(1, size(Hb1, 3));
    scale = 1e-5;
    for i = 1:size(Hb1, 3)
        Hbe = (squeeze(HbE1(:, :, i)) / scale);
        Hde = (squeeze(HdE1(:, :, i)) / scale);
        Rb = (squeeze(Hb1(:, :, i))'/ scale * squeeze(Hb1(:, :, i))/ scale);
        Rd = (squeeze(Hd1(:, :, i))'/ scale * squeeze(Hd1(:, :, i))/ scale);
        e = 1e-12;
        eta = 0.025 * sqrt(trace(Rb));
%         [EigvalueB, vectorB] = MaxEigenvector(Rb);
        EigvalueB = max(eig(Rb));
        gammaB = 0.5 * EigvalueB;
        gammaD = 0.05 * norm(Rd)/1e+2;
        alphaL = 1 / power((1 - eta / sqrt(EigvalueB)), 2);
        I = eye(size(Rd));
%         [EigvalueBD, vectorBD] = MaxEigenvector(pinv(Rd + gammaD * I) * Rb);
        EigvalueBD = max(eig(pinv(Rd + gammaD * I) * Rb));
%         w0 = ones(N, 1);
%         w0 = rand(N, 1, "like", 1i);
        w0 = rand(N, 1);
        alphaU = EigvalueBD * (w0') * (Rd + gammaD * I) * w0;
        
        % opt compute
        NoIterations = 20;
        alphaOpt = zeros(1, NoIterations);
        WOpt = zeros(N, N, NoIterations);
        for j = 1:NoIterations  
            % Calculate variables Dopt and αopt by solving (15).
            if j == 1
                alphaC = (alphaU - alphaL) .* rand(1,1) + alphaL;
%                 alphaC = (alphaU - alphaL) / 2;
            end
            cvx_begin 
                variable  W(N, N) complex 
                variable alphaopt 
                minimize(real(trace((Rd + gammaD * I) * W)))
                subject to
                    real(trace(Rb * W)) == alphaopt;
                    W == semidefinite(N, N);
                    real(alphaL) <= real(alphaopt);
                    real(alphaopt) <= real(alphaU);
                    real(power(eta, 2) * (trace(W)) + (sqrt(alphaC) - 1) + (alphaopt) * (1 / sqrt(alphaC) - 1)) <= 0;
            cvx_end
            WOpt(:, :, j) = W; 
            alphaOpt(j) = alphaopt;
            alphaC = alphaopt;
            if (j >= 2) && (trace((Rd + gammaD * I) * WOpt(:, :, j - 1)) - trace((Rd + gammaD * I) * WOpt(:, :, j)) <= e)
                break;  
            end
        end 
        wOpt = W(:, 1) / W(1, 1);
        PerformanceIndex = CaculateAC_NSDE;
        [AC(i)] = PerformanceIndex.AC(wOpt, Hbe, Hde, i);
    end
end
