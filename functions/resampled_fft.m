function[Hb, f, rirResampled] = resampled_fft(rir, nDFT, fsResampled, fs, directPath)
    rirResampled = 3*resample(rir, fsResampled, fs);
    if nargin > 4
        rirResampled((directPath + 5):end) = 0;
    end
    Y = fft(rirResampled, nDFT);
    P1 = Y(1:floor(nDFT/2) + 1);
    f = fsResampled*(0:floor(nDFT/2))/nDFT; 
    Hb= P1';
end