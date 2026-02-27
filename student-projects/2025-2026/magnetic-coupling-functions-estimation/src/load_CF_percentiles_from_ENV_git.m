function struct = load_CF_percentiles_from_ENV_git(folderPath,injectionType,building,excludeFolders,percValue)
    
    if nargin<4||isempty(excludeFolders);excludeFolders={};end
    excludeFolders=cellstr(excludeFolders);
    if nargin<5||isempty(percValue)||~isnumeric(percValue);error('loadCF:BadPercentile','percValue must be numeric');end
    percValue=percValue(:).';
    if ~isfolder(folderPath);error('loadCF:BadPath','Folder not found');end
    
    d=dir(folderPath);d=d([d.isdir]);d=d(~ismember({d.name},{'.','..'}));
    dateFolders={};
    for i=1:numel(d)
        if ~isempty(regexp(d(i).name,'^\d{8}$','once'));dateFolders{end+1}=d(i).name;end
    end
    if ~isempty(excludeFolders);dateFolders=setdiff(dateFolders,excludeFolders,'stable');end
    if isempty(dateFolders);error('loadCF:NoFolders','No valid folders');end
    dateFolders=sort(dateFolders);
    
    allFreq=[];allCF=[];
    
    for i=1:numel(dateFolders)
        currentFolder=fullfile(folderPath,dateFolders{i});
        filePattern=sprintf('*_%s_%s_CF.txt',injectionType,building);
        files=dir(fullfile(currentFolder,filePattern));
        for f=1:numel(files)
            filename=fullfile(currentFolder,files(f).name);
            opts=delimitedTextImportOptions("NumVariables",9);
            opts.DataLines=[2 Inf];opts.Delimiter="\t";
            opts.VariableNames=["FreqHz","CFmT","V3","V4","V5","V6","V7","V8","V9"];
            opts.SelectedVariableNames=["FreqHz","CFmT"];
            opts.VariableTypes=["double","double","string","string","string","string","string","string","string"];
            opts.ExtraColumnsRule="ignore";opts.EmptyLineRule="read";
            tbl=readtable(filename,opts);
            valid=tbl.FreqHz>0 & tbl.CFmT>0 & isfinite(tbl.CFmT);
            allFreq=[allFreq;tbl.FreqHz(valid)];
            allCF=[allCF;tbl.CFmT(valid)];
        end
    end
    
    if isempty(allFreq);error('loadCF:NoData','No valid CF data');end
    
    uFreq=unique(allFreq);
    nF=numel(uFreq);nP=numel(percValue);
    percCF=nan(nP,nF);
    
    for k=1:nF
        idx=allFreq==uFreq(k);
        for p=1:nP
            percCF(p,k)=prctile(allCF(idx),percValue(p));
        end
    end
    struct.freq=uFreq;
    struct.percCF=percCF;
    struct.allCF=allCF;
end