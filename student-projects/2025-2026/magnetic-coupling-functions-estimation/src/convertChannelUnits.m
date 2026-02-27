function out = convertChannelUnits(struct,convlist)
    out=struct;
    changes=false;
    for k=1:numel(convlist)
        token=strrep(convlist{k},' ','');
        parts=strsplit(token,'-->');
        if numel(parts)~=2
            error('convertChannelUnits:BadToken %s',token);
        end
        srcUnit=parts{1}; 
        dstUnit=parts{2};
        frBufs=fieldnames(struct.DATA);
        for bufIdx=1:numel(frBufs)
            currBufName=frBufs{bufIdx};
            for fn=fieldnames(struct.DATA.(currBufName).ChannelsUnitCurrent)'
                ch=fn{1};
                if strcmp(struct.DATA.(currBufName).ChannelsUnitCurrent.(ch),srcUnit)
                    fprintf("Converting channel '%s': ",ch)
                    switch [srcUnit '-->' dstUnit]
                        case 'nT-->T'
                            scale=1e-9;
                        case 'T-->nT'
                            scale=1e9;
                        case 'strain-->m'
                            scale=3000;
                        case 'm-->strain'
                            scale=1/3000;
                        otherwise
                            error('convertChannelUnits:UnsupportedConversion %s',token);
                    end
                    fprintf("'%s', ",token)
                    if ~strcmp(currBufName,'RAW') && isfield(out.DATA.(currBufName).Data,ch)
                        out.DATA.(currBufName).Data.(ch)=out.DATA.(currBufName).Data.(ch).*scale;
                        changes=true;
                    end
                    out.DATA.(currBufName).ChannelsUnitCurrent.(ch)=dstUnit;
                    out.DATA.(currBufName).ChannelsUnitHistory.(ch)=append(out.DATA.(currBufName).ChannelsUnitHistory.(ch),'-->',dstUnit);
                    fprintf("DONE.\n")
                end
            end
        end
    end
    if ~changes
        fprintf("No units were changed.\n")
    end
end