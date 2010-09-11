function model=molFilesToCDFfile(model,cdfFileName)
% Concatenates all the mol files in current folder into a cdf file
%
% Creates a cdf file, named  cdfFileNam, out of all the mol files in the 
% current folder. The cdf file can then be used with the web based
% implementation of the group contribution method to estimate the Standard
% Gibbs energy of formation for a batch of metabolite species
% The web-based implementation of this new group contribution method is
% available free at http://sparta.chem-eng.northwestern.edu/cgi-bin/GCM/WebGCM.cgi.
% The code checks for a .mol file with the filename prefix given by the
% abbreviation in the model, therefore, you should name your own mol files
% accordingly.
% 
%INPUT
% model.mets      cell array of metabolite abbreviations corresponding to the mol files
% cdfFileName     name of cdf file
%
%OUTPUT
% model.met(m).molFile  =1 when mol file found, otherwise =0
%
%WRITTEN OUTPUT
% cdfFileName.cdf    cdf file with all the mol files in order of the
%                   model metabolite abbreviations 
%
% Ronan M.T. Fleming

[nMet,nRxn]=size(model.S);
noMolFileCount = 0;

fid0=fopen([cdfFileName '.cdf'],'w');
for m=1:nMet
    metAbbr=model.mets{m};
    metAbbr=metAbbr(1:end-3);
    fid = fopen([metAbbr '.mol'],'r');
    if fid~=-1
        while 1
            tline = fgetl(fid);
            if ~ischar(tline)  
                break
            end
            fprintf(fid0,'%s\n',tline);
        end
        fclose(fid);
        fprintf(fid0,'%s\n','$$$$');
        model.met(m).molFile=1;
    else
        fprintf('%s\n',['No mol file for ' model.mets{m}]);
        model.met(m).molFile=0;
        noMolFileCount=noMolFileCount+1;
    end
end
fprintf('%s\n',['Percentage of metabolites without mol files: ' num2str(noMolFileCount/nMet)]);
fclose(fid0);