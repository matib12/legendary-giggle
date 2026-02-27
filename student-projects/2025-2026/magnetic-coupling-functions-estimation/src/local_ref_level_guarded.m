function ref = local_ref_level_guarded(x, idx0, Nguard, Nbin, method)
%LOCAL_REF_LEVEL_GUARDED  Local reference level excluding peak and guard bins
%
% ref = local_ref_level_guarded(x, idx0, Nguard, Nbin, method)
%
% INPUT:
%   x       : local spectrum (vector)
%   idx0    : central bin index (peak)
%   Nguard  : number of guard bins to exclude on each side
%   Nbin    : number of reference bins per side
%   method  : 'median' (default) or 'mean'
%
% OUTPUT:
%   ref     : reference level
%
% The function errors if there are not enough bins on both sides.

    % ------------------
    % Argument checks
    % ------------------
    if nargin < 5 || isempty(method)
        method = 'median';
    end

    % basic checks
    if ~isvector(x)
        error('x must be a vector.');
    end

    if ~isscalar(idx0) || idx0 <= 0 || idx0 ~= floor(idx0)
        error('idx0 must be a positive integer.');
    end

    if ~isscalar(Nguard) || Nguard < 0 || Nguard ~= floor(Nguard)
        error('Nguard must be a non-negative integer.');
    end

    if ~isscalar(Nbin) || Nbin <= 0 || Nbin ~= floor(Nbin)
        error('Nbin must be a positive integer.');
    end

    N = numel(x);

    if idx0 > N
        error('idx0 exceeds vector length.');
    end

    % ------------------
    % Check bin availability
    % ------------------
    minLeft  = idx0 - Nguard - Nbin;
    maxRight = idx0 + Nguard + Nbin;

    if minLeft < 1 || maxRight > N
        error(['Not enough bins around idx0. ', ...
               'Required: %d bins per side + %d guard bins.'], ...
               Nbin, Nguard);
    end

    % ------------------
    % Reference indices
    % ------------------
    idx_left  = (idx0 - Nguard - Nbin) : (idx0 - Nguard - 1);
    idx_right = (idx0 + Nguard + 1)    : (idx0 + Nguard + Nbin);

    idx_ref = [idx_left idx_right];
    if Nguard==0
        idx_ref = [idx_ref idx0];
        idx_ref = unique(idx_ref);
    end

    % ------------------
    % Reference level
    % ------------------
    switch lower(method)
        case 'median'
            ref = median(x(idx_ref), 'omitnan');
        case 'mean'
            ref = mean(x(idx_ref), 'omitnan');
        case 'max'
            ref = max(x(idx_ref),[], 'omitnan');
        otherwise
            error('Unknown method. Use ''median'' or ''mean'' or ''max''.');
    end

end