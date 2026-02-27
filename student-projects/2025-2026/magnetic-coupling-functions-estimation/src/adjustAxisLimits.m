function adjustAxisLimits(scaleFactor, axisType)
    if ~isscalar(scaleFactor) || scaleFactor < 0
        error('scaleFactor must be a non-negative scalar');
    end

    ax = gca;
    pl = findobj(ax,'-regexp','Type','line|scatter');
    
    
    % ------------------ axis scale ------------------
    switch axisType
        case 'X'
            scaleType = ax.XScale;
            dataAxis = 'XData';
        case 'Y'
            scaleType = ax.YScale;
            dataAxis = 'YData';
        otherwise
            error('axisType must be ''X'' or ''Y'' (case sensitive)');
    end

    % dmin = min(cellfun(@(x) min(x,[],'omitnan'), get(pl,dataAxis)));
    % dmax = max(cellfun(@(x) max(x,[],'omitnan'), get(pl,dataAxis)));
    dmin = min(arrayfun(@(h) min(h.(dataAxis),[],'omitnan'), pl));
    dmax = max(arrayfun(@(h) max(h.(dataAxis),[],'omitnan'), pl));

    % ------------------ linear scale ------------------
    if strcmp(scaleType,'linear')
        delta = dmax - dmin;

        if delta == 0
            delta = abs(dmin);
            if delta == 0
                delta = 1;
            end
        end

        newMin = dmin - delta * scaleFactor;
        newMax = dmax + delta * scaleFactor;

    % ------------------ logarithmic scale ------------------
    elseif strcmp(scaleType,'log')
        lmin=log10(dmin);
        lmax=log10(dmax);
        delta = lmax - lmin;

        if delta == 0
            delta = abs(lmin);
            if delta == 0
                delta = 1;
            end
        end

        newMin = 10^(lmin - delta * scaleFactor);
        newMax = 10^(lmax + delta * scaleFactor);

    else
        error('Unknown axis scale');
    end

    % ------------------ apply ------------------
    switch axisType
        case 'X'
            xlim([newMin newMax])
        case 'Y'
            ylim([newMin newMax])
    end

end