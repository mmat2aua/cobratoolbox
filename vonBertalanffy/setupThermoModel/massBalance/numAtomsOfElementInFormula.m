function N=numAtomsOfElementInFormula(formula,element)
% returns the number of atoms of a single element in a formula
%
% INPUT
% formula       formula in format Element then NumberOfAtoms with no spaces
% element       Abbreviation of element e.g. C or Mg
%
% OUTPUT
% N             number of atoms of this element in the formula provided
%
% Ronan Fleming 9 March 09
% Ronan Fleming 21 July 09 handles formulae with same first letter for
%                          element e.g.C & Co in 'C55H80CoN15O11'
% Ronan Fleming 18 Sept 09 handles composite formulae like C62H90N13O14P.C10H11N5O3.Co

if any(~isletter(element))
    disp(element)
    error('Element must not contain numbers')
end

if length(element)==1
    if ~strcmp(upper(element),element)
        disp(element)
        error('Single letter element must not be lower case')
    end
end

if ~isletter(formula(1))
    disp(formula)
    error('Formula format expected is element then number of elements, not the other way around')
end

zero='0';
indZero=strfind(formula,zero);
if ~isempty(indZero)
    if isletter(formula(indZero-1))
        % formula = strrep(formula, zero, 'O');
        error('Formula contains a zero with a letter preceeding it, represent oxygen with a the character O not zero')
    end
end

%exceptions due to metabolites with wierd formula (Human reconstruction)
formula = strrep(formula, 'FULLRCO', '');
formula = strrep(formula, 'FULLR2CO2', '');
formula = strrep(formula, 'FULLRFULLR2CO', '');
formula = strrep(formula, 'C3H52FULLR3CO2', '');
formula = strrep(formula, 'C20H29OFULLR2CO', '');
formula = strrep(formula, 'XCO2', '');

if ischar(element) && ischar(formula)
    %handle composite formulae like C62H90N13O14P.C10H11N5O3.Co
    indDots=findstr('.',formula);
    if ~isempty(indDots)
        N=0;
        lastDot=1;
        for z=1:length(indDots)+1
            if z <= length(indDots)
                formulaD=formula(lastDot:indDots(z)-1);
                lastDot=indDots(z);
            else
                formulaD=formula(lastDot:end);
            end
            
            %now inner iteration to count elements
            start=strfind(formulaD,element);
            %be sure to distinguish between two elements with the same starting letter
            %in a formulaD, e.g. C & Co
            xlt=length(start);
            flt=length(formulaD);
            if xlt>1
                x=xlt;
                while x>0
                    if start(x) + 1 <= flt
                        if isletter(formulaD(start(x)+1))
                            if strcmp(formulaD(start(x)+1),lower(formulaD(start(x)+1)))
                                x=x-1;
                            end
                        else
                            %letter with number afterward
                            start=start(x);
                            break;
                        end
                    else
                        %single letter with no number afterward taken as correct
                        %element e.g. something strange like ...Co10C
                        start=start(x);
                        break;
                    end
                end
            end
            if ~isempty(start)
                %if this is a single letter element
                if start~=length(formulaD) && length(element)==1
                    %and if next character is a letter
                    if isletter(formulaD(start+1))
                        %and if that letter is lower case
                        if strcmp(formulaD(start+1),lower(formulaD(start+1)))
                            %then this is not a match e.g. C vs Ca
                            start=[];
                        end
                    end
                end
            end
            if ~isempty(start)
                %initialise
                ND=0;
                %start looking for number after end of element character
                start=start+length(element);
                if start>length(formulaD)
                    %only one atom
                    ND=1;
                else
                    %molecule
                    p=0;
                    %read until no numbers left
                    while (start+p)<=length(formulaD)
                        %stop if this is this a not a number
                        if isempty(str2num(formulaD(start+p)))
                            break;
                        else
                            %add another character that is a number
                            p=p+1;
                        end
                    end
                    if p==0
                        %entirely text formulas, e.g. HCN, with a single atom of
                        %particular element
                        ND=1;
                    end
                end
                if ND~=1
                    %molecule nummer
                    ND=str2num(formulaD(start:start+p-1));
                end
            else
                %no atoms of this element in this formulaD
                ND=0;
            end
            %count the elements up from different parts of the composite
            %formula
            N=N+ND;
        end
    else
        start=strfind(formula,element);

        %be sure to distinguish between two elements with the same starting letter
        %in a formula, e.g. C & Co
        xlt=length(start);
        flt=length(formula);
        if xlt>1
            x=xlt;
            while x>0
                xlast=x;
                if start(x) + 1 <= flt
                    if isletter(formula(start(x)+1))
                        if strcmp(formula(start(x)+1),lower(formula(start(x)+1)))
                            x=x-1;
                        end
                    else
                        %letter with number afterward
                        start=start(x);
                        break;
                    end
                else
                    %single letter with no number afterward taken as correct
                    %element e.g. something strange like ...Co10C
                    start=start(x);
                    break;
                end
                if xlast==x
                    fprintf('%s\n',formula);
                    error('Stuck on formula.')
                end
            end
        end
        if ~isempty(start)
            %if this is a single letter element
            if start~=length(formula) && length(element)==1
                %and if next character is a letter
                if isletter(formula(start+1))
                    %and if that letter is lower case
                    if strcmp(formula(start+1),lower(formula(start+1)))
                        %then this is not a match e.g. C vs Ca
                        start=[];
                    end
                end
            end
        end
        if ~isempty(start)
            %initialise
            N=0;
            %start looking for number after end of element character
            start=start+length(element);
            if start>length(formula)
                %only one atom
                N=1;
            else
                %molecule
                p=0;
                %read until no numbers left
                while (start+p)<=length(formula)
                    %stop if this is this a not a number
                    if isempty(str2num(formula(start+p)))
                        break;
                    else
                        %add another character that is a number
                        p=p+1;
                    end
                end
                if p==0
                    %entirely text formulas, e.g. HCN, with a single atom of
                    %particular element
                    N=1;
                end
            end
            if N~=1
                %molecule nummer
                N=str2num(formula(start:start+p-1));
            end
        else
            %no atoms of this element in this formula
            N=0;
        end

    end
else
    if ~ischar(element)
        disp(element)
        error('Element must be given by a variable of class char')
    else
        disp(element)
        error('Formula must be given by a variable of class char')
    end
end