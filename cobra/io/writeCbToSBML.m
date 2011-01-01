function writeCbToSBML(model,fileName)
%writeCbToSBML write a COBRA model to SBML
%
% writeCbToSBML(model,fileName)
%
%INPUTS
% model         COBRA model structure
% fileName      Name of xml file output
%
        sbmlModel = cobra_struct_to_sbml_struct(model);
        OutputSBML(sbmlModel, fileName);
