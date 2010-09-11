function GetMolFileFromKeggID(Name,OutputFileName)
% function connects to Kegg database and retrieves the mol files for each compound with a Kegg ID
%
% Name              model structure or xls file name
% OutputFileName    Filename for output
%
% Ines Thiele June 29, 2010
%

if isstruct(Name)
    model = Name;
    Mode = 1;
else
    strExcelFile = Name;
    Mode = 2;
end

warning off;

if Mode == 2
% The KEGG database 
strLinks='KeggID'; % 

fprintf('Reading URLs from spreadsheet\n')
[numeric,txt,raw]=xlsread(strExcelFile);
[nrows,ncols]=size(raw);
iCol = strmatch(strLinks,regexprep(txt(1,:),' ',''),'exact');

elseif Mode == 1
    raw{:,1} = model.metKeggID; 
    iCol = 1;
end

i=2; % Starting row

fprintf('Checking URLs ...\n')
cnt = 1;
while i<nrows
    % Parse data from Excel table
    fprintf('%d/%d\n', i, nrows)
    if ~isempty(raw{i,iCol}) &&  length(char(raw{i,iCol}))>1
        strURL = strcat('http://www.genome.jp/dbget-bin/www_bget?-f+m+compound+',raw{i,iCol});
        n=0;
        s=urlread(strURL);
        if ~isempty(s)
            FoundMolFile{cnt,1} = raw{i,iCol};
            dlmwrite(OutputFileName,strcat('New Compound: ',raw{i,iCol}),'-append','delimiter','');
            dlmwrite(OutputFileName,s,'-append','delimiter','');
            cnt = cnt + 1;
        end
    end
    i=i+1;
end

