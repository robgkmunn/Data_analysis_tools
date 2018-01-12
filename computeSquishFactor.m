function [rho_max, lambda] = computeSquishFactor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script will run the expansion_factor script on a spreadsheet of data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all windows and clean things up
close all; clc; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spreadsheet = uigetfile('*.xlsx','Choose input file');
k = 1;
% get data ft
[~,~,data] = xlsread(spreadsheet,k);
spreadsheet_short = strsplit(spreadsheet,'.'); spreadsheet_short = spreadsheet_short{1};
% load the session, tetrode, and unit to load the right data
row1 = data(1,:);
Session_base = data(2:end,(strcmp(row1, 'Sessions') == 1));
Session_squish = data(2:end,strcmp(row1, 'Sessions_1') == 1);
Tetrode_base = data(2:end,find(strcmp(row1, 'Tetrode') == 1,1,'first'));
Tetrode_squish = data(2:end,find(strcmp(row1, 'Tetrode_1') == 1,1,'first'));
Unit_base = data(2:end,(strcmp(row1, 'Unit') == 1));
Unit_squish = data(2:end,(strcmp(row1, 'Unit_1') == 1));
Box = data(2:end,(strcmp(row1, 'Box') == 1));

%initialize vectors to hold the rho_max & lambda
rho_max = nan(length(Session_base),1);
lambda = nan(length(Session_base),1);

% file1 = Spike file, open field
% file2 = Position file, open field
% file3 = Spike file, squish
% file4 = Position file, squish

for n = 1:length(Session_base)
    boxSize = Box{n};
    file2 = strcat(Session_base{n},'_pos.mat');
    file1 = strcat(Session_base{n},'_T',num2str(Tetrode_base{n}),'C',num2str(Unit_base{n}),'.mat');
    file4 = strcat(Session_squish{n},'_pos.mat');
    file3 = strcat(Session_squish{n},'_T',num2str(Tetrode_squish{n}),'C',num2str(Unit_squish{n}),'.mat');
    [rho_max(n), yShift(n), xShift(n) lambda(n)] = squish_factor_translation_RM_v3(file1, file2, file3, file4, boxSize);
end
keyboard
