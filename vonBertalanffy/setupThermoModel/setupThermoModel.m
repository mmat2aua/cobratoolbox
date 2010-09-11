function [modelT,directions,solutionThermoRecon]=setupThermoModel(model,metAbbrAlbertyAbbr,...
    metGroupCont,Alberty2006,temp,PHA,IS,CHI,biomassRxnAbbr,...
    rxnAbbrDatabaseID,metBoundsFile,rxnBoundsFile,Legendre,useKeqData,...
    nStdDevGroupCont,cumNormProbCutoff,figures,printToFile)
% Thermodynamically constrains a COBRA model.
%
% This function takes as imput an ordinary COBRA model (such as from sbml.xml)
% then uses thermodynamic data to setup a thermodynamic COBRA model.
%
%INPUT
% model
% model.biomassRxnAbbr      abbreviation of biomass reaction
%
% Alberty2006  Basic data on the metabolite species that make
%              up a reactant, compiled by Robert A. Alberty,
%              Massachusetts Institute of Technology.
%              In Print: Robert A. Alberty, Biochemical Thermodynamics: 
%              Applications of Mathematica. John Wiley & Sons, 2006. p391-395
%              Online: BasicBioChemData3.nb
%              http://library.wolfram.com/infocenter/MathSource/5704/
%              Only species of interest in the range pH 5 to 9 are included.
%
% Alberty2006 is a structure with two fields for each metabolite:
% Alberty2006.abbreviation(i) Alberty reactant abbreviation
% Alberty2006.basicData(i) cell array with 4 columns: dGf0,dHf0,charge,#Hydrogens 
%
% metAbbrAlbertyAbbr    mapping from model metabolite primary key to
%                       primary key of reactants in Alberty2006
%
% metGroupCont  Group contribution data, derived by submitting metabolite 
%               .mol files to the web based implementation, as described in:
%               M. D. Jankowski, C. S. Henry, L. J. Broadbelt, and V. Hatzimanikatis.
%               Group contribution method for thermodynamic analysis of complex metabolic networks. 
%               Biophys J, 95(3):1487-1499, Aug 2008.
% 
% metGroupCont has the following fields
% metGroupCont(m).abbreviation                      metabolite abbreviation
% metGroupCont(m).delta_G_formation                 estimate of standard Gibbs energy of formation    
% metGroupCont(m).delta_G_formation_uncertainty     standard error in estimate of standard Gibbs energy of formation
% metGroupCont(m).formulaMarvin                     metabolite formula of molfile (Marvin)
% metGroupCont(m).chargeMarvin                      metabolite charge of molfile (Marvin)
% metGroupCont(m).groupContribution_pH              glass electrode pH of metabolite molfile
% metGroupCont(m).groupContribution_file            fileName with origin of group contribution data
%
% temp              temperature 273.15 K to 313.15 K
%
% PHA.*             glass electrode pH, between 5-9, in compartment defined by letter *
%                   *   compartment
%                   c   cytoplasm
%                   p   periplasm
%                   e   extracellular environment
%                   m   mitochondria
%                   n   nucleus
%                   r   endoplasmic reticulum
%                   g   Golgi apparatus
%                   l   lysosome
%                   x   peroxisome
%
% IS.*              ionic strength (0 - 0.35M) in compartment defined by letter *
% 
% CHI.*             electrical potential (mV) in compartment defined byletter *
%
% rxnAbbrDatabaseID   n x 2 cell array of reaction abbreviation and corresponding database ID 
%
%
% OPTIONAL INPUT
%
% metBoundsFile     filename with upper & lower bounds on metabolite 
%                   concentrations (Molar)
%                   See 'model_met_bounds.txt' for format required
% rxnBoundsFile     filename with upper & lower bounds on fluxes 
%                   (mmol gDw-1 hr-1)
%                   See 'model_rxn_bounds.txt' for format required
%
% Legendre              {(1),0} Legendre Transformation of dGf0?
% useKeqData            {(1),0} Use dGf0 back calculated from Keq?
% nStdDevGroupCont      {1} number of standard deviations of group contribution
%                       uncertainty to be used for each metabolite
% cumNormProbCutoff     {0.1} positive real number between 0 and 0.5 that
%                       specifies to tolerance when there is uncertainty in group
%                       contribution estimates.
% figures               {0}, 1 = create figures
% printToFile           {0}, 1 = print out to log files
%
%OUTPUT
% modelT        a thermodynamic constraints based model
%
% modelT.gasConstant       Gas Constant
% modelT.faradayConstant   Faraday Constant
% modelT.temp              Temp
%
% modelT.met(m).dGft0                standard transformed Gibbs energy of formation(kJ/mol)
% modelT.met(m).dGft0Source          origin of data, Keq or groupContFileName.txt 
% modelT.met(m).dGft0Keq             standard transformed Gibbs energy of formation(kJ/mol) 
% modelT.met(m).dGft0GroupCont            group. cont. estimate of standard transformed Gibbs energy of formation(kJ/mol)
% modelT.met(m).dGf0GroupCont             group. cont. estimate of standard Gibbs energy of formation(kJ/mol)
% modelT.met(m).dGf0GroupContUncertainty  group. cont. uncertainty in estimate of standard Gibbs energy of formation (kJ/mol) 
% modelT.met(m).mf           mole fraction of each species within a pseudoisomer group
% modelT.met(m).aveZi        average charge
% modelT.met(m).chi          electrical potential
% modelT.met(m).aveHbound    average number of protons bound to a reactant
%
% modelT.lb_reconThermo      lower bound from amalgamation of reconstruction
%                            and thermodynamic assignment directions
% modelT.ub_reconThermo      upper bound from amalgamation of reconstruction
%                            and thermodynamic assignment directions
%
% directions    a structue of boolean vectors with different directionality 
%               assignments where some vectors contain subsets of others
%
% qualitatively assigned directions 
%       directions.fwdReconBool
%       directions.revReconBool
%       directions.reversibleReconBool
%
% qualitatively assigned directions using thermo in preference to
% qualitative assignments but using qualitative assignments where
% thermodynamic data is lacking
%       directions.fwdReconThermoBool
%       directions.revReconThermoBool
%       directions.reversibleReconThermoBool
%
% reactions that are qualitatively assigned by thermodynamics
%       directions.fwdThermoOnlyBool
%       directions.revThermoOnlyBool
%       directions.reversibleThermoOnlyBool
%
% qualtiative -> quantiative changed reaction directions 
%       directions.ChangeReversibleFwd
%       directions.ChangeReversibleRev
%       directions.ChangeForwardReverse
%       directions.ChangeForwardReversible
%
% subsets of forward qualtiative -> reversible quantiative change
%       directions.ChangeForwardReversibleBool_dGfGC
%       directions.ChangeForwardReversibleBool_dGfGC_byConcLHS
%       directions.ChangeForwardReversibleBool_dGfGC_byConcRHS
%       directions.ChangeForwardReversibleBool_dGfGC_bydGt0
%       directions.ChangeForwardReversibleBool_dGfGC_bydGt0LHS
%       directions.ChangeForwardReversibleBool_dGfGC_bydGt0Mid
%       directions.ChangeForwardReversibleBool_dGfGC_bydGt0RHS
%       directions.ChangeForwardReversibleBool_dGfGC_byConc_No_dGt0ErrorLHS
%       directions.ChangeForwardReversibleBool_dGfGC_byConc_No_dGt0ErrorRHS
%
%       directions.cumNormProbCutoff
%       directions.ChangeForwardForwardBool_dGfGC
%
% solutionThermoRecon       FBA solution with thermodynamic in preference to 
%                           reconstruction directions 
%                           (see setThermoReactionDirectionality.m)
%
% DEPENDENCIES ON OTHER COBRA TOOLBOX FUNCTIONS
% findRxnIDs.m
% printRxnFormula.m
%
% OPTIONAL DEPENDENCIES ON OTHER PACKAGES (i.e. more accurate)
% Multiple Precision Toolbox for MATLAB by Ben Barrowes (mptoolbox_1.1)
% http://www.mathworks.com/matlabcentral/fileexchange/6446
%
% COBRA TOOLBOX solvers e.g. solveCobraLP.m and solveCobraLPCPLEX.m for
% conflict resolution for an infeasible LP.
%
% This code is part of the vonBertalanffy package for Systems Biothermodynamics.
% Ludwig von Bertalanffy
% http://en.wikipedia.org/wiki/Ludwig_von_Bertalanffy
%
% Ronan M.T. Fleming, Sept 2010.

if ~exist('Legendre','var')
    Legendre=1;
end
if ~exist('useKeqData','var')
    useKeqData=1;
end
if ~exist('nStdDevGroupCont','var')
    nStdDevGroupCont=1;
end
if ~exist('cumNormProbCutoff','var')
    cumNormProbCutoff=0.1;
end
if ~exist('figures','var')
    figures=0;
end
if figures==1
    close all
end
if ~exist('printToFile','var')
    printToFile=0;
end

%all possible compartments
numChar=1;
[allMetCompartments,uniqueCompartments]=getCompartment(model.mets,numChar);

if ~isfield(model,'biomassRxnAbbr')
    if ~exist('biomassRxnAbbr')
        fprintf('\n%s\n','...checkObjective');
        objectiveAbbr=checkObjective(model);
        fprintf('%s\n',['Asumming objective is ' objectiveAbbr]);
        model.biomassRxnAbbr=objectiveAbbr;
    else
        model.biomassRxnAbbr=biomassRxnAbbr;
    end
end

model=findSExRxnInd(model);
fprintf('\n');

if printToFile
    %create a folder to contain the directionality report
    if ~isfield(model,'description')
        model.description='';
    end
    folderName=[model.description '_thermoDirectionality'];
    mkdir(folderName)
    fullFolderName=[pwd '/' folderName];
    cd(fullFolderName);
end

[nMet,nRxn]=size(model.S);
if printToFile==1
    [imBalancedMass,imBalancedCharge,imBalancedBool] = checkMassChargeBalance(model,[],-1);
else
    [imBalancedMass,imBalancedCharge,imBalancedBool] = checkMassChargeBalance(model,[],1);
end
if isempty(imBalancedBool)
    imBalancedBool=false(nRxn,1);
end

if 0
    if ~isempty(imBalancedMass)
        error('Internal reactions are not mass balanced.')
    end
    if ~isempty(imBalancedCharge)
        error('Internal reactions are not charge balanced.')
    end
else
    if ~isempty(imBalancedMass)
        warning('Internal reactions are not mass balanced.')
    end
    if ~isempty(imBalancedCharge)
        warning('Internal reactions are not charge balanced.')
    end
end

%Physico-Chemical Constants (Energies are expressed in kJ/mol)
gasConstant = 8.314472/1000; % kJ K-1 mol-1
faradayConstant = 96.48;     % Coloumb mol-1
model.gasConstant=gasConstant;
model.faradayConstant=faradayConstant;
%stamp model with intensive thermodynamic constants
model.PHA=PHA;
model.IS=IS;
model.temp=temp;

%optionally map reaction database ID's to model structure
if exist('rxnAbbrDatabaseID','var')
    model=mapDatabaseID(model,rxnAbbrDatabaseID);
end

fprintf('\n%s\n','...convert to COBRA v2 style format');
%DASHES NOT changed to UNDERSCORES in abbreviations
model=convertToCobraV2(model);

%adjust for reactions involving co2 and succinate dehydrogenase reaction
model=thermodynamicAdjustmentToStoichiometry(model);
[nMet,nRxn]=size(model.S);

%readjust the length of the boolean vector indicating unbalanced reactions
imBalancedBool2=false(nRxn,1);
imBalancedBool2(1:length(imBalancedBool),1)=imBalancedBool;
imBalancedBool=imBalancedBool2;


% Sets any reactions that have equal upper bounds to be forward reaction
fprintf('%s\n','Checking for reactions that have equal upper and lower bounds...');
for n=1:nRxn
    if model.lb(n)==0 && model.ub(n)==0
%         error(['Reaction ' model.rxns{n} ' lb=ub=zero'])
        warning(['Reaction ' model.rxns{n} ' lb=ub=zero']);
        fprintf('%s\n',['Reaction ' model.rxns{n} ' ' model.rxn(n).officialName ' set to forward']);
        model.lb(n)=0;
        model.ub(n)=1000;
        model.rxn(n).directionality='forward';
        model.rxn(n).regulationStatus='Off';
    else
        model.rxn(n).regulationStatus='On';
    end
end
    
%new reactions may have been added so update model size
[nMet,nRxn]=size(model.S);

%add constraint sense, including for new rows, if any.
if ~isfield(model,'csense')
    model.csense(1:nMet)='E';
end

% load Ecoli_metName_albertyAbbreviation;
nAlbertyMet=size(metAbbrAlbertyAbbr,1);
%replace any underscores in metabolite abbreviations with dashes
for m=1:length(metGroupCont)
    x = strfind(metGroupCont(m).abbreviation,'_');
    if x~=0
        metGroupCont(m).abbreviation(x)='-';
    end
end

fprintf('%s\n','...replace any underscores in metabolite abbreviations with dashes');
for m=1:nAlbertyMet
    abbr=metAbbrAlbertyAbbr{m,2};
    x = strfind(abbr,'_');
    if x~=0
        abbr(x)='-';
        metAbbrAlbertyAbbr{m,2}=abbr;
    end
end

fprintf('\n%s\n','...mapping Alberty abbreviations to model using metAbbrAlbertyAbbr.');
if printToFile
    fid=fopen('metabolites_with_no_Alberty_abbreviation.txt','w');
end
for m=1:nMet
    got=0;
    for a=1:nAlbertyMet
        metAbbr=model.mets{m};
        %dont include compartment
        metAbbr=metAbbr(1:end-3);
        if strcmp(metAbbr,metAbbrAlbertyAbbr{a,2})
            if ~isempty(metAbbrAlbertyAbbr{a,3})
                model.met(m).albertyAbbreviation=metAbbrAlbertyAbbr{a,3};
                got=1;
            end
            break;
        end
    end
    if got==0
        metAbbr=model.mets{m};
        if strcmp(metAbbr(end-2:end),'[c]')
            if printToFile==0
                fprintf('%s\t%s\t%20s\t%s\n','No Alberty abbreviation for metabolite', int2str(m),model.mets{m},model.met(m).officialName);
            else
                fprintf(fid,'%s\t%s\t%s\t%s\n','No Alberty abbreviation for metabolite', int2str(m),model.mets{m},model.met(m).officialName);
            end
        end
    end
end
if printToFile
    fclose(fid);
end

fprintf('\n%s\n','...assignment of Group Contribution data to model using metGroupCont structure.');
%assign Group Contribution data to model using metGroupCont
 fprintf('%s\n','Energies are expressed in kJ mol^-1')
nMetGroupCont=length(metGroupCont);
NaNdGf0GCMetBool=false(nMet,1);
metGroupContAbbr=cell(nMetGroupCont,1);
for a=1:nMetGroupCont
    metGroupContAbbr{a,1}=metGroupCont(a).abbreviation;
end
%check if compartment identifier is also in metabolite abbreviation
if strcmp(metGroupCont(1).abbreviation(end),']')
    %check if metabolite abbreviation in model matches any in group
    %contribution data
    for m=1:nMet
        bool=strcmp(model.mets{m},metGroupContAbbr);
        if ~any(bool)
            %mark as missing
            NaNdGf0GCMetBool(m,1)=1;
            model.met(m).dGf0GroupCont=NaN;
            model.met(m).dGf0GroupContUncertainty=NaN;
            model.met(m).formulaMarvin=NaN;
            model.met(m).chargeMarvin=NaN;
            model.met(m).groupContribution_pH=NaN;
            model.met(m).groupContribution_file=NaN;
        else
            if nnz(bool)>1
                error([metAbbr ': duplicated abbreviation[*] in group contribution data']);
            else
                %chemical standard chemical potential is  independent of compartment
                %TODO - group contribution estimate specific to pH
                %Energies are expressed in kJ mol^-1 so convert from  kCal
                model.met(m).dGf0GroupCont=metGroupCont(bool).delta_G_formation*(8.314472/1.987);
                model.met(m).dGf0GroupContUncertainty=metGroupCont(bool).delta_G_formation_uncertainty*(8.314472/1.987);
                model.met(m).formulaMarvin=metGroupCont(bool).formulaMarvin;
                model.met(m).chargeMarvin=metGroupCont(bool).chargeMarvin;
                model.met(m).groupContribution_pH=metGroupCont(bool).pH;
                model.met(m).groupContribution_file=metGroupCont(bool).file;
            end
        end
    end
else
    for m=1:nMet
        if strcmp(model.met(m).abbreviation,'damval[c]');
            pause(eps)
        end
        
        %check if metabolite abbreviation in model matches any in group
        %contribution data
        metAbbr=model.met(m).abbreviation;
        metAbbr=metAbbr(1:end-3);
        bool=strcmp(metAbbr,metGroupContAbbr);
        if ~any(bool)
            %mark as missing
            NaNdGf0GCMetBool(m,1)=1;
            
            model.met(m).dGf0GroupCont=NaN;
            model.met(m).dGf0GroupContUncertainty=NaN;
            model.met(m).formulaMarvin=NaN;
            model.met(m).chargeMarvin=NaN;
            model.met(m).groupContribution_pH=NaN;
            model.met(m).groupContribution_file=NaN;
        else
            bool=find(bool);
            bool=bool(1);
            if nnz(bool)>1
                error([metAbbr ': duplicated abbreviation in group contribution data']);
            else
                %chemical standard chemical potential is  independent of compartment
                %TODO - Assign pH to input for OpenBabel to get InChi strings specific
                %to pH
                %Energies are expressed in kJ mol^-1 so convert from  kCal
                model.met(m).dGf0GroupCont=metGroupCont(bool).delta_G_formation*(8.314472/1.987);
                model.met(m).dGf0GroupContUncertainty=metGroupCont(bool).delta_G_formation_uncertainty*(8.314472/1.987);
                model.met(m).formulaMarvin=metGroupCont(bool).formulaMarvin;
                model.met(m).chargeMarvin=metGroupCont(bool).chargeMarvin;
                model.met(m).groupContribution_pH=metGroupCont(bool).pH;
                model.met(m).groupContribution_file=metGroupCont(bool).file;
            end
        end
    end
end


% Apparent glass electrode pH is not the same as real pH for thermodynamic calculations.
% Given the experimental glass electrode measurement of pH, this function returns
% the real pH to be used for thermodynamic calculations, pHr = -log10[H+], 
% by subtracting the effect of the ion atmosphere around H+ which 
% reduces its activity coefficient below unity.
% See p49 Alberty 2003
fprintf('\n%s\n','...realpH');
for p=1:length(uniqueCompartments)
    if isfield(PHA,uniqueCompartments{p,1})
        %if comparing this program against Albertys tables then use an apparent pH, i.e. pHa,
        %with a thermodynamic pH, i.e. pHr, that is equivalent
        compareAgainstAlbertysTables=1; %changed
        if compareAgainstAlbertysTables
            [pHr,pHAdjustment]=realpH(PHA.(uniqueCompartments{p,1}),temp,IS.(uniqueCompartments{p,1}));
            PHA.(uniqueCompartments{p,1})=PHA.(uniqueCompartments{p,1})+pHAdjustment;
        else
            fprintf('\n%s\n','If comparing this data with Albertys tables, note that they use real pH.');
        end
        
        [pHr,pHAdjustment]=realpH(PHA.(uniqueCompartments{p,1}),temp,IS.(uniqueCompartments{p,1}));
        %real thermodynamic pH
        PHR.(uniqueCompartments{p,1})=pHr;
    end
end

%Incorportate the electrochemical potential across the membrane by adding a
%component to the standard chemical potential of external protons:
%"Escherichia coli Glutamate- and Arginine-Dependent Acid Resistance Systems Increase Internal pH and
%Reverse Transmembrane Potential" by Hope Richard and John W. Foster*
%Data from Table 1
%Given a medium at pH 7 the exterior of the cell has an electrical
%potential 90mV lower than the exterior. Interior pH of 7.8. Culture temp
%of 310.15K.

%By default, make a Legendre transformation for electrical potential
LegendreCHI=1;

fprintf('\n%s\n','...assignThermoToModel');
%assign Alberty data to model, or Legendre transformed Henry data if no data 
%on metabolite is available from Alberty
%this function is where most of the detailed physical chemistry implemented
%by Alberty in mathematica, is implemented in matlab.
model=assignThermoToModel(model,Alberty2006,temp,PHR,IS,CHI,uniqueCompartments,NaNdGf0GCMetBool,Legendre,LegendreCHI,useKeqData,printToFile);

fprintf('\n%s\n','...setCommonZeroStandardGibbsEnergyOfFormation');
% Set all the exceptional metabolites to have a common baseline.
% i.e. Standard transformed Gibbs energies of reactants with the baseline adjusted
% for certain paired cofactors e.g. fad & fadh2, such that the
% difference between the two is the same as in Albertys data but
% the absolute values are consistent with the group contribution data 
model=setCommonZeroStandardGibbsEnergyOfFormation(model);

%balance the protons in each reaction given the number of Hydrogens bound
%to each reactant calculated thermodynamically using assignThermoToModel.m
%note that where a reaction involves two compartments, the products are
%assumed to dissociate. This needs more investigation in order to
%definitively prove that this assumption does not result in net production
%of protons due to complex shuttling of metabolites over and back across a
%membrane. Follow up work on this is planned. - Ronan
% e.g. iAF1260
% pH balancing of protons:
% ThermoRecon directions (Relaxed all conflicting directions). Growth:  0.755689
% c.f. growth with reconstruction directions:                           0.722621
% Without pH balancing of protons:
% ThermoRecon directions (Relaxed all conflicting directions). Growth:  0.775954
% c.f. growth with reconstruction directions:                           0.736701

if 0
fprintf('\n%s\n','...pHbalanceProtons');
model=pHbalanceProtons(model,imBalancedBool);
end

%plot statistics on the number of reactants with significant
%non-predominant mole fractions >0.05
if figures
    moleFractionStats(model)
end

fprintf('\n%s\n','... readMetRxnBoundsFiles');
%assign bounds on metabolite concentrations and fluxes

% setDefaultConc            sets default bounds on conc [1e-5,0.02]     
setDefaultConc=1;
% setDefaultFlux            sets all reactions reversible [-1000,1000]
setDefaultFlux=0;% set to zero since we use the bounds given by

%locations of the files with the bounds on metabolite concentration 
%and reaction flux
if ~exist('metBoundsFile','var')
%     metBoundsFile='model_met_bounds.txt';
%    metBoundsFile='Schuetz_met_data.txt';
    metBoundsFile='Bennet_Glucose_Aerobic.txt';
%     metBoundsFile='Bennet_Glucose_Aerobic_Cofactor.txt';
%    metBoundsFile='Bennet_Glucose_Aerobic_Cofactor_ATPs.txt';
%      metBoundsFile='Bennet_Glucose_Aerobic_Cofactor_NAD.txt';
end
if ~exist('rxnBoundsFile','var')
%     rxnBoundsFile='model_rxn_bounds.txt';
    rxnBoundsFile='Ecoli_glucose_aerobic_rxn_bounds.txt';
end
% metBoundsFile=[];
% rxnBoundsFile=[];
%assign bounds to model +/- read in data from flat file 
model=readMetRxnBoundsFiles(model,setDefaultConc,setDefaultFlux,metBoundsFile,rxnBoundsFile);

%Special concentration bounds for o2, co2, h2o and h
for m=1:nMet
    abbr=model.mets{m};
    abbrShort=abbr(1:end-3);
    compartment=abbr(end-1);
    
    %water concentration assumed to be one molar p107 Alberty 2003
    if strcmp('h2o',abbrShort)
        model.met(m).concMin=0.99;%55.5062-1;;
        model.met(m).concMax=1;%55.5062+1;;
    end

    if strcmp('co2',abbrShort)
        if strcmp('co2[c]',abbr)
            model.met(m).concMin=10e-8;
            model.met(m).concMax=0.0014;
        else
            model.met(m).concMin=0.0001;
            model.met(m).concMax=0.0001;
        end
    end
    
    %From Henry et al
    % oxygen concentration selected for the media is 8.2E-6 M
    % and the oxygen concentration in the cell cannot exceed the concentration
    % in the media, the bounds on the oxygen concentration in the cell
    % were set from 10E-7 M to 8.2E-6M
    if strcmp('o2',abbrShort)
        if strcmp('o2[c]',abbr)
            model.met(m).concMin=10e-8;
            model.met(m).concMax=8.2e-6;
        else
            model.met(m).concMin=8.2e-8;
            model.met(m).concMax=8.2e-6;
        end
    end
    %hydrogen ion concentration uses real pH (not apparent pH)
    if strcmp('h',abbrShort)
        model.met(m).concMin=10^-PHR.(compartment);
        model.met(m).concMax=10^-PHR.(compartment);
    end
end

fprintf('\n%s\n','... deltaG0concFluxConstraintBounds');
%Set reaction directionality bounds from thermodynamic data:
%set up bounds  on metabolite chemical potential
%use to set upper & lower thermodynamic bounds on internal fluxes
%******************first pass***************
model=deltaG0concFluxConstraintBounds(model,Legendre,LegendreCHI,figures,nStdDevGroupCont);

fprintf('\n%s\n','...standardGibbsFormationEnergyStats');
% figures=0;
[nKeq,nGC,nNone]=standardGibbsFormationEnergyStats(model,figures);
% figures=1;

fprintf('\n%s\n','...directionalityCheckAll');
%check the thermodynamically feasible directions with respect to the
%reconstruction directions
if printToFile
    %print out problematic reactions in summary table format for paper
    printToTable=1;
    % cumNormProbCutoff     {0.1} positive real number between 0 and 0.5 that
    %                       specifies to tolerance when there is uncertainty in group
    %                       contribution estimates.
else
    printToTable=0;
end

%create new standard Gibbs Free Energy based on the geometric mean of each
%metabolites concentration range
thorStandard=0; %zero for first pass is adivisable since growth may be infeasible  
directions=directionalityCheckAll(model,cumNormProbCutoff,thorStandard,printToFile,printToTable,figures);
%boolean vectors indexing the reaction directionalities according to
%different criteia
model.directions=directions;

%test the functionality of the model if cobra toolbox installed
if exist('solveCobraLP','file')==0
    fprintf('\n')
    fprintf('%s\n','No LP solver configured with a COBRA toolbox installation.');
    fprintf('%s\n','See http://gcrg.ucsd.edu/Downloads/Cobra_Toolbox');
    fprintf('%s\n','The second pass assignment of reaction directionality requires an LP solver.');
    solutionRecon=[];
    solutionThermoRecon=[];
    model1=[];
else
    global CBTLPSOLVER
    %FBA with qualitatively assigned directions using thermo in preference to
    % qualitative assignments but using qualitative assignments where
    % thermodynamic data is lacking
    % model.lb_reconThermo              lower bounds from dGtMin/dGtMax and recon
    %                                   directions if thermo data missing
    % model.ub_reconThermo              upper bounds from dGtMin/dGtMax and recon
    %                                   directions if thermo data missing
    modelD=model;
    modelD.lb=model.lb_reconThermo;
    modelD.ub=model.ub_reconThermo;
    if strcmp(CBTLPSOLVER,'cplex_direct')
        printLevel=1;
        basisReuse=0;
        %set conflict resolve to 1 if you want a log file to be printed which
        %will help to identify the source of the infeasibility. The code in
        %cplex uses Chinnek's heuristics which dont always work but often they
        %are useful.
        %http://www.sce.carleton.ca/faculty/chinneck/chinneck_pub.shtml
        conflictResolve=0;
        if conflictResolve==1
            %set the lower bound on the objective to a minimum value to force
            %biomass production
            modelD.lb(strcmp(modelD.rxns,modelD.biomassRxnAbbr))=0.1;
        end
        contFunctName=[];
        %minimising the norm of a flux vector will eliminate flux around loops
        %but this is only guarunteed to give a thermodynamically feasible flux
        %if there are no upper or lower bounds on internal reactions.
        minNorm=1e-6;
        [solutionThermoRecon,modelD]=solveCobraLPCPLEX(modelD,printLevel,basisReuse,conflictResolve,contFunctName,minNorm);
    else
        printLevel=3;
        primalOnlyFlag=[];
        solver=CBTLPSOLVER;
        modelD.A=modelD.S; modelD.osense = -1;
        changeOK = changeCobraSolverParams('LP','printLevel',printLevel);
        solutionThermoRecon = solveCobraLP(modelD);
        modelD=rmfield(modelD,'A');model=rmfield(modelD,'osense');
    end
    clear modelD
    fprintf('\n%s\t%g\n','Growth with amalgamation of recon & thermo directions: ',solutionThermoRecon.obj);
    %******************second pass starts here ***************
    secondPass=0;
    if secondPass
        %requires manual curation
        [model,solutionThermoRecon,solutionRecon,model1]=secondPassDirectionalityAssignment(model);
    end
end

fprintf('\n%s\n','...readableCobraModel');
%make the model readable
model=readableCobraModel(model);

if Legendre==1 && LegendreCHI
    fprintf('\n%s\n','N.B. Thermodynamic properties calculated using a Legendre transform for pH and electrical potential.');
else
    fprintf('\n%s\n','N.B. Thermodynamic properties have not been calculated using a Legendre transform.');
end

%move out of folder
cd ..
if printToFile
    fprintf('\n%s\n',['Directionality report in folder: ' folderName]);
end
%change name of model
modelT=model;