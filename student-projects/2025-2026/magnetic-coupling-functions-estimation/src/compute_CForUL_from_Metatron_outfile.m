function out = compute_CForUL_from_Metatron_outfile(destinationStruct,filename,FFTlength,overlapPerc,SNRthrHrec,SNRthrSensor,percentiles)
    [~, fname, ~] = fileparts(filename);
    tokens = regexp(fname,'^(?<type>[^_]+)_.*_(?<building>[A-Z]{3})-\d+$','names');
    building = tokens.building;
    type     = tokens.type;

    HrecChannelName     = 'Hrec_hoft_2_200Hz';
    HrecRawChannelName  = 'Hrec_hoft_20000Hz';

    switch type
        case 'MagneticLine'
            switch building
                case 'NEB'
                    channelBuildingSubName = 'NE';
                    asdSquareModulusName = "ENV_NE_MAG";
                case 'WEB'
                    channelBuildingSubName = 'WE';
                    asdSquareModulusName = "ENV_WE_MAG";
                case 'CEB'
                    channelBuildingSubName = 'CEB';
                    asdSquareModulusName = "ENV_CEB_MAG";
                otherwise
                    error('compute_CForUL_from_Metatron_outfile:BadBuilding', "Passed a not valid building: '%s'", building)
            end
            sensorChannelNames = {['ENV_' channelBuildingSubName '_MAG_N'],...
                ['ENV_' channelBuildingSubName '_MAG_V'],...
                ['ENV_' channelBuildingSubName '_MAG_W']};
        otherwise
            error('compute_CForUL_from_Metatron_outfile:BadAnalysisType', "Passed a not valid analysis type: '%s'", type)
    end
    if ~size(fieldnames(destinationStruct),1)==0
        fprintf("Data loading for channels: ")
        allNamesForPrint=[sensorChannelNames,HrecChannelName,HrecRawChannelName,destinationStruct.CHANNEL];
        for s = 1:numel(allNamesForPrint)
            fprintf("'%s', ",allNamesForPrint{s})
        end
        fprintf("skipped, data already present.\n")
        out=destinationStruct;
    else
        out = importMetatronMagneticLinesInjections(filename);
        NoiseChannelName = ['ENV_' out.CHANNEL];
        out = importDataMetratronMagneticInjections(out,'rds',[sensorChannelNames,HrecChannelName,NoiseChannelName]);
    end
    out = convertChannelUnits(out,{'nT-->T', 'strain-->m'});
    out = estimateSpectrums(out,'RDS',[sensorChannelNames,{HrecChannelName}],FFTlength,overlapPerc);
    % out = estimateSpectrums(out,'RAW',sensorChannelNames,FFTlength,overlapPerc);
    out = buildQuadratureASD(out,'RDS',sensorChannelNames,asdSquareModulusName);
    % % % out = buildQuadratureASD_RAW(out,sensorChannelNames,asdSquareModulusName);
    out = computeCouplingFunctionOrUpperLimit(out,'RDS','all',HrecChannelName,asdSquareModulusName,1,2,15,'max',SNRthrHrec,SNRthrSensor);
    out=computeNoiseProjections(out,out.InjParams.GPS_Start(1)-1200,600, sensorChannelNames,100,50,percentiles);
end

