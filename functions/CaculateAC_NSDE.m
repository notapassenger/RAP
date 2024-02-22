function PerformanceIndex = CaculateAC_NSDE
    PerformanceIndex.AC_NSDE = @AC_NSDE;
    PerformanceIndex.AC = @AC;
end
function[AC, NSDE] =  AC_NSDE(w, pd, Hbe, Hde, i)
    % input: w, pd, HbE, HdE
    Mb = size(Hbe, 1);
    Md = size(Hde, 1);
    % compute AC
    mu = ((Md*(w')*(Hbe')*Hbe*w) / (Mb*(w')*(Hde')*Hde*w));
    AC = real(10*log10(mu));
    % compute NSDE
    PM = (power(norm(Hbe*w-pd), 2) / power(norm(pd), 2));
    NSDE = 10*log10(PM);
end

function[AC] = AC(w, Hbe, Hde, i)
    Mb = size(Hbe, 1);
    Md = size(Hde, 1);
    % compute AC
    mu = ((Md*(w')*(Hbe')*Hbe*w) / (Mb*(w')*(Hde')*Hde*w));
    AC = real(10*log10(mu));
end