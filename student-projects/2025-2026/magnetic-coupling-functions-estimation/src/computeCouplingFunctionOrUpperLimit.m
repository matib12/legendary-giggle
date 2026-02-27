function struct = computeCouplingFunctionOrUpperLimit(struct,frBufName,freqsel,hrec,sensor,NbinInj,NbinHrecQuiet,NbinSensorOffsetQuiet,method,SNRthrHrec,SNRthrSensor)
    if ~isfield(struct.DATA.(frBufName),'ASD')||~isstruct(struct.DATA.(frBufName).ASD),error('computeCouplingFunction:MissingASD','ASD missing or not a struct');end
    if ~isfield(struct.DATA.(frBufName),'ASDf'),error('computeCouplingFunction:MissingASDf','ASDf missing');end
    if ~isfield(struct.DATA.(frBufName).ASD,hrec),error('computeCouplingFunction:MissingHrec','Hrec %s missing in ASD',hrec);end
    if ~isfield(struct.DATA.(frBufName).ASD,sensor),error('computeCouplingFunction:MissingSensor','Sensor %s missing in ASD',sensor);end
    if NbinInj<0, error('computeCouplingFunction:NbinInvalid','Number of bin to consider for frequency selection lower than zero');end
    if ~any(strcmp(method,["median","mean","max"])), error('computeCouplingFunction:InvalidMethod',"Method '%s' not valid, use 'median' or 'mean' or 'max'",method);end

    Nh  = size(struct.DATA.(frBufName).ASD.(hrec),1);
    Ns  = size(struct.DATA.(frBufName).ASD.(sensor),1);
    
    if ischar(freqsel)
        if ~strcmp(freqsel,'all') || Nh~=Ns
            error('computeCouplingFunction:BadFreqSel','Invalid freq selection');
        end
        freqsel=1:length(struct.InjParams.Frequency);
    else
        if any(freqsel<=0) || any(freqsel>Nh) || any(freqsel>Ns)
            error('computeCouplingFunction:BadFreqIdx','Requested rows do not exist');
        end
    end

    fprintf("Estimating CF between '%s' and '%s', SNRthrHrec = %.2f, SNRthrSensor = %.2f.\nComputed frequencies: ", hrec,sensor,SNRthrHrec,SNRthrSensor)
    
    if isfield(struct.DATA.(frBufName),'CForUL'),struct.DATA.(frBufName)=rmfield(struct.DATA.(frBufName),'CForUL');end
    struct.DATA.(frBufName).CForUL.SNRthrHrec=SNRthrHrec;
    struct.DATA.(frBufName).CForUL.SNRthrSensor=SNRthrSensor;

    for fridx=freqsel
        frInj=struct.InjParams.Frequency(fridx);
        [~, idxFrInj] = min(abs(struct.DATA.(frBufName).ASDf - frInj));
        fprintf("%.3f Hz, ",frInj)

        currHrec        = struct.DATA.(frBufName).ASD.(hrec)(fridx,:);
        currSensor      = struct.DATA.(frBufName).ASD.(sensor)(fridx,:);

        CF_Hrec_Inj     = local_ref_level_guarded(currHrec,    idxFrInj,0,NbinInj,method);
        CF_Hrec_Quiet   = local_ref_level_guarded(currHrec,    idxFrInj,NbinInj,NbinHrecQuiet,method);
        CF_Sensor_Inj   = local_ref_level_guarded(currSensor,  idxFrInj,0,NbinInj,method);
        CF_Sensor_Quiet = local_ref_level_guarded(currSensor,  idxFrInj,NbinInj+NbinSensorOffsetQuiet,NbinHrecQuiet*4,method);

        SNR_Hrec   = CF_Hrec_Inj   / CF_Hrec_Quiet;
        SNR_Sensor = CF_Sensor_Inj / CF_Sensor_Quiet;

        if ~isfield(struct.DATA.(frBufName),'CForUL') || ~isfield(struct.DATA.(frBufName).CForUL,'SearchedFrequencies')
            struct.DATA.(frBufName).CForUL.SearchedFrequencies=frInj;
        else
            struct.DATA.(frBufName).CForUL.SearchedFrequencies=[struct.DATA.(frBufName).CForUL.SearchedFrequencies frInj];
        end

        if SNR_Hrec > SNRthrHrec && SNR_Sensor > SNRthrSensor
            CF = sqrt((CF_Hrec_Inj^2-CF_Hrec_Quiet^2)/(CF_Sensor_Inj^2-CF_Sensor_Quiet^2));
            if ~isfield(struct.DATA.(frBufName).CForUL,'CF')
                struct.DATA.(frBufName).CForUL.CF=CF;
                struct.DATA.(frBufName).CForUL.freqCF=frInj;
                struct.DATA.(frBufName).CForUL.SNR_CF_Hrec=SNR_Hrec;
                struct.DATA.(frBufName).CForUL.SNR_CF_Sensor=SNR_Sensor;
            else
                struct.DATA.(frBufName).CForUL.CF=[struct.DATA.(frBufName).CForUL.CF CF];
                struct.DATA.(frBufName).CForUL.freqCF=[struct.DATA.(frBufName).CForUL.freqCF frInj];
                struct.DATA.(frBufName).CForUL.SNR_CF_Hrec=[struct.DATA.(frBufName).CForUL.SNR_CF_Hrec SNR_Hrec];
                struct.DATA.(frBufName).CForUL.SNR_CF_Sensor=[struct.DATA.(frBufName).CForUL.SNR_CF_Sensor SNR_Sensor];
            end
        else
            UL = CF_Hrec_Quiet / sqrt(CF_Sensor_Inj^2-CF_Sensor_Quiet^2);
            if ~isfield(struct.DATA.(frBufName).CForUL,'UL')
                struct.DATA.(frBufName).CForUL.UL=UL;
                struct.DATA.(frBufName).CForUL.freqUL=frInj;
                struct.DATA.(frBufName).CForUL.SNR_UL_Hrec=SNR_Hrec;
                struct.DATA.(frBufName).CForUL.SNR_UL_Sensor=SNR_Sensor;
            else
                struct.DATA.(frBufName).CForUL.UL=[struct.DATA.(frBufName).CForUL.UL UL];
                struct.DATA.(frBufName).CForUL.freqUL=[struct.DATA.(frBufName).CForUL.freqUL frInj];
                struct.DATA.(frBufName).CForUL.SNR_UL_Hrec=[struct.DATA.(frBufName).CForUL.SNR_UL_Hrec SNR_Hrec];
                struct.DATA.(frBufName).CForUL.SNR_UL_Sensor=[struct.DATA.(frBufName).CForUL.SNR_UL_Sensor SNR_Sensor];
            end
        end
        % close all
        % reset(groot)
        % figure('Position',[100 100 1200 800])
        % subplot(1,2,1)
        % hold on
        % loglog(struct.DATA.(frBufName).ASDf(idxFrInj-1:idxFrInj+1), currHrec(idxFrInj-1:idxFrInj+1),'or','MarkerFaceColor', 'red', 'DisplayName','INJ bins')
        % loglog(struct.DATA.(frBufName).ASDf([idxFrInj-1-NbinHrecQuiet:idxFrInj-2, idxFrInj+2:idxFrInj+1+NbinHrecQuiet]), currHrec([idxFrInj-1-NbinHrecQuiet:idxFrInj-2, idxFrInj+2:idxFrInj+1+NbinHrecQuiet]),'og', 'MarkerFaceColor', 'green', 'DisplayName','QUIET bins')
        % loglog(struct.DATA.(frBufName).ASDf, currHrec,'.k','MarkerFaceColor','k','DisplayName','Inj')
        % loglog(struct.DATA.(frBufName).ASDf,repelem(CF_Hrec_Inj,length(struct.DATA.(frBufName).ASDf)),'--', 'DisplayName','CF Hrec Inj')  
        % loglog(struct.DATA.(frBufName).ASDf,repelem(CF_Hrec_Quiet,length(struct.DATA.(frBufName).ASDf)),'--', 'DisplayName','CF Hrec Quiet')  
        % xlim([frInj-0.1 frInj+0.1])
        % legend
        % grid minor
        % set(gca(),'XScale','log','Yscale','log')
        % 
        % subplot(1,2,2)
        % hold on
        % loglog(struct.DATA.(frBufName).ASDf(idxFrInj-1:idxFrInj+1), currSensor(idxFrInj-1:idxFrInj+1),'or','MarkerFaceColor', 'red', 'DisplayName','INJ bins')
        % loglog(struct.DATA.(frBufName).ASDf([idxFrInj-1-NbinHrecQuiet*4-15:idxFrInj-2-15, idxFrInj+2+15:idxFrInj+1+NbinHrecQuiet*4+15]), currSensor([idxFrInj-1-NbinHrecQuiet*4-15:idxFrInj-2-15, idxFrInj+2+15:idxFrInj+1+NbinHrecQuiet*4+15]),'og', 'MarkerFaceColor', 'green', 'DisplayName','QUIET bins')
        % loglog(struct.DATA.(frBufName).ASDf, currSensor,'.k','MarkerFaceColor','k','DisplayName','Inj')
        % loglog(struct.DATA.(frBufName).ASDf,repelem(CF_Sensor_Inj,length(struct.DATA.(frBufName).ASDf)),'--', 'DisplayName','CF Sensor Inj')  
        % loglog(struct.DATA.(frBufName).ASDf,repelem(CF_Sensor_Quiet,length(struct.DATA.(frBufName).ASDf)),'--', 'DisplayName','CF Sensor Quiet')  
        % 
        % xlim([frInj-0.1 frInj+0.1])
        % legend
        % grid minor
        % set(gca(),'XScale','log','Yscale','log')
        % sgtitle(sprintf("SNR Hrec = %.2f, SNR Sensor = %.2f",SNR_Hrec,SNR_Sensor))
        % % a=frInj;
        % % strtmp=sprintf("tmpImages/%s_frInj_%dHz%d.png",sensor,floor(a),floor(mod(a,floor(a))*100));
        % % exportgraphics(gcf, strtmp, 'Resolution', 600);
    end
    fprintf(' Done.\n')
end