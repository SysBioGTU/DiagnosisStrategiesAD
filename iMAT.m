%%%%%% iMAT algorithm to reconstruct personalized metabolic models %%%%%%%% 
clear all; clc;
% Human-GEM (Version 1.11.0) was retrieved from the GitHub page
% (https://github.com/SysBioChalmers/Human-GEM) and saved as .mat file.
load('cobraHuman.mat')

model = simplifyModel(cobraHuman); % Firstly remove boundary metabolites
 
model.lb(7719)= -1000;  % lower bound for glucose uptake
model.ub(7719)= -0.01;  % upper bound for glucose uptake
model.lb(7733)= -1000;  % lower bound for oxygen uptake
model.ub(7733)= -0.01;  % upper bound for oxygen uptake
model.lb(13350)=0.0001; % lower bound for humna biomass reaction
 
model.A= model.S; % Model_T.S; % Stochiomatrix format for gurobi showen as A instead of S
model.A=sparse(model.A);
model.rhs=model.b;
model.obj=model.c; % objective function
model.sense=['=']; 
model.vtype='C' ;

% Normalized and covariate adjusted transcriptome data (from Neff et al.)
% was read here.
[num,txt,row] = xlsread('file_name.xlsx');

expressionData.value = num;
expressionData.gene = txt(2:end,1);  

% to decide thresholds
[num1,txt1,row1] = xlsread('all(AD_NCI).xlsx');
exp_data1= num;
y1 = quantile(mean(exp_data1') , 7);

%initCobraToolbox
% changeCobraSolver('gurobi', 'all'); 
% addpath('tests')
% addpath('utils')
% addpath('wrappers')
% addpath('plotting')


N=size(expressionData.value, 2);
Result_ALL = zeros(length(model.rxns) , N);
mymodel_Backup_0 = model;

for i=1:N
    Temp_Model = mymodel_Backup_0;
    Temp_expressionData.gene = expressionData.gene;
    Temp_expressionData.value = expressionData.value( : , i);
    [Temp_expressionRxns_AD Temp_parsedGPR] = mapExpressionToReactions(Temp_Model, Temp_expressionData);
    Temp_new_model = iMAT(Temp_Model, Temp_expressionRxns_AD, y1(2) , y1(6));
    clear Temp_Model Temp_expressionData Temp_parsedGPR Temp_expressionRxns_AD;

    k = ismember(mymodel_Backup_0.rxns, Temp_new_model.rxns);
    Result_ALL=[Result_ALL k];
    clear k;
    clear Temp_new_model; 
end

xlswrite('binary_models.xlsx',Result_ALL);

save('analysis.mat')
