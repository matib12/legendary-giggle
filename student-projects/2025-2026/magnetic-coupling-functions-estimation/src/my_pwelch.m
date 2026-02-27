function [Pxx, f] = my_pwelch(x, fs, npoints, overlapPerc, outType)
%MY_PWELCH Welch PSD/ASD with fixed window and user-defined overlap
%
%   [Pxx, f] = my_pwelch(x, fs, npoints, overlapPerc)
%   [Pxx, f] = my_pwelch(x, fs, npoints, overlapPerc, outType)
%
%   INPUT:
%     x           : data vector (1D)
%     fs          : sampling frequency [Hz]
%     npoints     : segment length and nfft
%     overlapPerc : overlap percentage (0 <= overlapPerc < 100)
%     outType     : (optional) 'asd' -> returns sqrt(PSD)
%
%   OUTPUT:
%     Pxx         : PSD [unit^2/Hz] or ASD [unit/sqrt(Hz)]
%     f           : frequency vector [Hz]

    arguments
        x (:,1) double
        fs (1,1) double {mustBePositive}
        npoints (1,1) {mustBePositive, mustBeInteger}
        overlapPerc double = []
        outType (1,:) char = 'psd'
    end

    % Force column vector
    x = x(:);

    % Welch parameters
    window   = hanning(npoints);
    noverlap = round(npoints * overlapPerc / 100);
    nfft     = npoints;

    % Compute PSD using Welch method
    [Pxx, f] = pwelch(x, window, noverlap, nfft, fs);

    % Output format selection
    if strcmpi(outType, 'asd')
        Pxx = sqrt(Pxx);
    end

end