function output = importMetatronMagneticLinesInjections(filename)
    fid = fopen(filename,'r');
    assert(fid>0, 'Cannot open file');

    fgetl(fid); %To read first line and go to the one with values
    headerValueLine  = fgetl(fid);
    fclose(fid);

    headerNames = {'Quiet_GPS_Start','Quiet_Duration','INJECTION_TYPE','CHANNEL'};
    headerValues = strsplit(strtrim(headerValueLine));

    output=struct();

    for k = 1:numel(headerNames)
        num = str2double(headerValues{k});
        if ~isnan(num)
            output.(headerNames{k}) = num;
        else
            output.(headerNames{k}) = headerValues{k};
        end
    end

    opts = detectImportOptions(filename);
    opts.DataLines = [4 Inf];
    opts.VariableNames = {'GPS_Start','Duration','Frequency','Amplitude'};
    data = readtable(filename, opts);

    output.InjParams = data;
end