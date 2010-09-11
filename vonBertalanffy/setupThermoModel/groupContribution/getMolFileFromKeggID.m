function molFiles = getMolFileFromKeggID(model)
% function connects to Kegg database and retrieves the mol files for each compound with a Kegg ID
%
%model.metKEGGID


fprintf('Checking URLs ...\n')
cnt = 1;
for i=1:size(model.S,1)
        strURL = strcat('http://www.genome.jp/dbget-bin/www_bget?-f+m+compound+',model.metKEGGID{i});
        s=urlread(strURL);
        if ~isempty(s)
            molFiles{i,1} = s;
            %dlmwrite(OutputFileName,strcat('New Compound: ',s),'-append','delimiter','');
            dlmwrite([model.mets{i}(1:end-3) '.mol'],s,'-append','delimiter','');
        end
end

