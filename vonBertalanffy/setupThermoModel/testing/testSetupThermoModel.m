%script for testing functionality of vonBertalanffy COBRA toolbox extension

%setup the path to the stripts and data
initVonBertalanffy

%modelToTest='iCore';
%modelToTest='iAF1260';
modelToTest='ReconX';

switch modelToTest
    case 'iAF1260'
        
        %name of biomass reaction
        biomassRxnAbbr='Ec_biomass_iAF1260_core_59p81M';
        
        load([vonBdir filesep 'setupThermoModel' filesep 'experimentalData' filesep 'reconstructions' filesep 'iAF1260_flux1.mat']);
        %group contrbuion data for E. coli
        load([vonBdir filesep 'setupThermoModel' filesep 'experimentalData' filesep 'groupContribution' filesep 'metGroupCont_Ecoli_iAF1260.mat']);
        
        model.S(952,350)=1;%one reaction needing mass balanced in iAF1260
        
        % default thermodynamic parameters
        if 0
            temp=310.15;
            PHA.c=7.7;
            PHA.p=7.7;
            PHA.e=7.7;
            IS.c=0.25;
            IS.p=0.25;
            IS.e=0.25;
            CHI.c=0;
            CHI.p=90; %milliVolts
            CHI.e=90;
        else
            temp=310.15;
            PHA.c=7.7;
            PHA.p=5;
            PHA.e=5;
            IS.c=0.25;
            IS.p=0.15;
            IS.e=0.15;
            CHI.c=0;
            CHI.p=90; %milliVolts
            CHI.e=90;
        end
        
    case 'iCore'
        %Ecoli core
        biomassRxnAbbr='Biomass_Ecoli_core_w/GAM';
        %Ecoli core
        load([vonBdir filesep 'setupThermoModel' filesep 'experimentalData' filesep 'reconstructions' filesep 'ecoli_core_xls2model.mat']);
        %group contrbuion data for E. coli
        load([vonBdir filesep 'setupThermoModel' filesep 'experimentalData' filesep 'groupContribution' filesep 'metGroupCont_Ecoli_iAF1260.mat']);
        
        
        % default thermodynamic parameters
        if 0
            temp=310.15;
            PHA.c=7.7;
            PHA.p=7.7;
            PHA.e=7.7;
            IS.c=0.25;
            IS.p=0.25;
            IS.e=0.25;
            CHI.c=0;
            CHI.p=90; %milliVolts
            CHI.e=90;
        else
            temp=310.15;
            PHA.c=7.7;
            PHA.p=5;
            PHA.e=5;
            IS.c=0.25;
            IS.p=0.15;
            IS.e=0.15;
            CHI.c=0;
            CHI.p=90; %milliVolts
            CHI.e=90;
        end
    case 'ReconX'
        %biomass reaction
        biomassRxnAbbr='biomass_reaction';
        if 1
            load([vonBdir filesep 'setupThermoModel' filesep 'experimentalData' filesep 'reconstructions' filesep 'ReconX.mat']);
            if isfield(model,'metCharge')
                model.metCharges=model.metCharge;
                rmfield(model,'metCharge')
            end
        else
            load('/home/common/IEM/IEMdirectionality/data/ReconX/Hulda/ReconX')
            
            model=newModel;
        end
        model.description='ReconX';
        %group contrbuion data for E. coli
        load([vonBdir filesep 'setupThermoModel' filesep 'experimentalData' filesep 'groupContribution' filesep 'metGroupCont_ReconX.mat']);
        
        %thermodynamic parameters
        temp=310.15;
        
        %log10(hydrogen ion activity)
        PHA.c=7.35;
        PHA.e=7;
        PHA.g=6.5;
        PHA.l=4.8;
        PHA.m=8;
        PHA.n=7.35;
        PHA.r=7.07;
        PHA.x=7;
        
        %ionic strength
        someIonicStrength=0.1;
        IS.c=someIonicStrength;
        IS.e=someIonicStrength;
        IS.g=someIonicStrength;
        IS.l=someIonicStrength;
        IS.m=someIonicStrength;
        IS.n=someIonicStrength;
        IS.r=someIonicStrength;
        IS.x=someIonicStrength;
        
        %milliVolts
        CHI.c=0;
        CHI.e=0;
        CHI.g=0;
        CHI.l=0;
        CHI.m=150; %Wan et al.
        CHI.n=0;
        CHI.r=0;
        CHI.x=0;
        
        
end

pause(eps)

if 1
    %parameters for routine setup WITH printouts and figures
    metBoundsFile=[];
    rxnBoundsFile=[];
    Legendre=1;
    useKeqData=1;
    nStdDevGroupCont=1;
    cumNormProbCutoff=0.2;
    figures=1;
    printToFile=1;
else
    %other choice of parameters for routine setup without printouts or
    %figures
    metBoundsFile=[];
    rxnBoundsFile=[];
    Legendre=1;
    useKeqData=1;
    nStdDevGroupCont=1;
    cumNormProbCutoff=0.2;
    figures=0;
    printToFile=0;
end

rxnAbbrDatabaseID=[];

%setup a thermodynamically constrained model of metabolism
[modelT,directions,solutionThermoRecon]=setupThermoModel(model,metAbbrAlbertyAbbr,...
    metGroupCont,Alberty2006,temp,PHA,IS,CHI,biomassRxnAbbr,...
    rxnAbbrDatabaseID,metBoundsFile,rxnBoundsFile,Legendre,useKeqData,...
    nStdDevGroupCont,cumNormProbCutoff,figures,printToFile);