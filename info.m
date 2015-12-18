%% Sequential blinks removal


%%
% date = '20151212';
% cellDescrip = 'RyR_Cy3B'; %'RyR_Cy3B' 'Cav12_AF647' RyR_AF555 FO_AF647
% channel = 'gree';
% runNum = 2; %1 3 5 7 9 11
%%
date = ['20151217'];
cellDescrip = 'Cav12_AF647'; %'RyR_Cy3B' 'Cav12_AF647' RyR_AF555 FO_AF647
channel = 'red';
runNum = [1]; %1 3 5 7 9 11
%%
% date = '20151118';
% cellDescrip = 'IgM_MZ_Cy3B'; %'RyR_Cy3B' 'Cav12_AF647' RyR_AF555 FO_AF647
% channel = 'green';
% runNum = [1 2 3 10]; %1 3 5 7 9 11 %1 3 5 7 9 11
%%
% date = ['20151116'; '20151117'];
% cellDescrip = ['Cav_Cy3B';'RyR_Cy3B';'RyR_Cy3B']; %'RyR_Cy3B' 'Cav12_Alexa647' RyR_AF555 FO_AF647
% channel = 'green';
% runNum = zeros(3,6);
% runNum(1,1:3) = 1:3;
% runNum(2,1:5) = 1:5; %1 3 5 7 9 11
% runNum(3,1:6) = 1:6;

%%

% date = ['20151030'; '20151105';'20151106'];
% cellDescrip = ['RyR_AF647';'Cav_AF647';'Cav_AF647']; %'RyR_Cy3B' 'Cav12_Alexa647' RyR_AF555 FO_AF647
% channel = 'red';
% runNum = zeros(3,6);
% runNum(1,1:3) = 1:3;
% runNum(2,1:5) = 1:5; %1 3 5 7 9 11
% runNum(3,1:6) = 1:6;
%%
% date = ['20151116'; '20151117'];
% cellDescrip = ['IgD_FO_AF647';'IgD_MZ_AF647']; %'RyR_Cy3B' 'Cav12_Alexa647' RyR_AF555 FO_AF647
% channel = 'red';
% runNum = zeros(2,4);
% runNum(1,1:4) = 1:4;
% runNum(2,1:4) = 1:4; %1 3 5 7 9 11
%%
% date = ['20151116'; '20151117'];
% cellDescrip = ['CD45_FO_Cy3B';'CD45_MZ_Cy3B']; %'RyR_Cy3B' 'Cav12_Alexa647' RyR_AF555 FO_AF647
% channel = 'green';
% runNum = zeros(2,4);
% runNum(1,1:4) = 1:4;
% runNum(2,1:4) = 1:4; %1 3 5 7 9 11
%%
% date = ['20150923'; '20150924';'20151016'];
% cellDescrip = ['IgM_FO_AF647';'IgM_MZ_AF647';'IgM_FO_AF647']; %'RyR_Cy3B' 'Cav12_Alexa647' RyR_AF555 FO_AF647
% channel = 'red';
% runNum = zeros(3,6);
% runNum(1,1:5) = 1:5;
% runNum(2,1:6) = 1:6; %1 3 5 7 9 11
% runNum(3,1:4) = 1:4;
%% Example
%{
date =      ['20140703';        '20140702';     '20140630';         '20140630';       '20140612';     '20140610'  ];
cellDescrip = [ 'Cav12-Phos ';    'Cav12      ';    'RyR        ';     'RyR-Phos   ';  'Microtubule' ; 'Microtubule' ];
runNum = zeros(6,21);
runNum(1, 1:8) = 1:8;
runNum(2, 1:11) = 1:11;
runNum(3, 1:15) = 1:15;
runNum(4, 16:21) = 16:21;
runNum(5, 1:10) = 1:10;
runNum(6, 1:10) = 1:10;
%}
