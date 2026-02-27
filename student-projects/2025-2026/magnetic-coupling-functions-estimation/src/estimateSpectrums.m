function struct = estimateSpectrums(struct,frBufName,chanlist,FFTlength,overlapPerc)
    if ~strcmp(frBufName,'RAW')
        if ~isfield(struct.DATA.(frBufName), 'Data') || ~isstruct(struct.DATA.(frBufName).Data)
            error('estimateSpectrums:MissingData\nField struct.Data is missing or not a struct.');
        end
    end
    for j = 1:numel(chanlist)
        ch = chanlist{j};
        if ~strcmp(frBufName,'RAW')
            if ~isfield(struct.DATA.(frBufName).Data, ch)
                error('estimateSpectrums:MissingChannel\nChannel "%s" not found in struct.Data.', ch);
            end
            if isempty(struct.DATA.(frBufName).Data.(ch))
                error('estimateSpectrums:EmptyChannel\nChannel "%s" in struct.Data is empty.', ch);
            end
        end
    end

    fs=struct.DATA.(frBufName).fs;
    npoints=FFTlength*fs;
    
    fprintf("Computing ASDs for '%s': ",frBufName)
    for j=1:numel(chanlist)
        ch=chanlist{j};
        fprintf("'%s', ", ch)
        if ~strcmp(frBufName,'RAW')
            for i=1:size(struct.InjParams,1)
                x=struct.DATA.(frBufName).Data.(ch)(i,:); 
                x=x(isfinite(x));
                [P,f]=my_pwelch(x,fs,npoints,overlapPerc,'asd');
                if j==1 && i==1
                    struct.DATA.(frBufName).ASDf=f';
                end
                if i==1
                    struct.DATA.(frBufName).ASD.(ch)=nan(size(struct.InjParams,1),numel(P));
                end
                struct.DATA.(frBufName).ASD.(ch)(i,:)=P;            
            end
        end
        struct.DATA.(frBufName).ASDUnits.(ch)=append(struct.DATA.(frBufName).ChannelsUnitCurrent.(ch),'/sqrt{Hz}');
    end
    fprintf(" DONE.\n")
end
