function out = importDataMetratronMagneticInjections(out, framebuffer, channelsName)
    switch framebuffer
        case {'rds', 'rds.ffl', '/virgoData/ffl/rds.ffl'}
            framebuffer = '/virgoData/ffl/rds.ffl';
            fsample     = 50;
            chansuffix  = '_50Hz';
            frBufName   = 'RDS';
        case {'raw', 'raw.ffl', '/virgoData/ffl/raw.ffl'}
            framebuffer = '/virgoData/ffl/raw.ffl';
            fsample     = 20000;
            chansuffix  = '';
            frBufName   = 'RAW';
        otherwise
            error('importDataMetratronMagneticInjections:InvalidFrameBuffer\nFramebuffer "%s" not supported or not implemented yet.', framebuffer);
    end

    datastruct.DataFrameBuffer  = framebuffer;
    datastruct.DataChannelsName = channelsName;
    datastruct.fs               = fsample;
    
    for j = 1:numel(channelsName)
        datastruct.Data.(channelsName{j}) = nan(size(out.InjParams,1),max(out.InjParams.Duration)*fsample); %Create nan matrix (Nrow x relativeDuration*fsample) for INJECTION DATA
    end

    fprintf("Getting channels from '%s': ",framebuffer)

    for j = 1:numel(channelsName)
        chName = channelsName{j};
        fprintf("'%s'", chName)
        for k = 1:size(out.InjParams,1)
            [data,units] = getChannelWithChannelUnit(framebuffer, [chName chansuffix], out.InjParams{k,1}, out.InjParams{k,2});
            datastruct.Data.(chName)(k,1:length(data)) = data(:);
            if k==1
                datastruct.ChannelsUnitOriginal.(chName)=units;
            end
        end
        fprintf(", ")
    end
    [datastruct.ChannelsUnitCurrent,datastruct.ChannelsUnitHistory]=deal(datastruct.ChannelsUnitOriginal);
    out.DATA.(frBufName)=datastruct;
    fprintf("DONE.\n")
end