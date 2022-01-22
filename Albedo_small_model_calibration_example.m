clc;clear all;close all;

%notes:
%Noah-MP SNEQV0 may differ from that in this small model, because this
%SNEQV0 value is not output, rather calculated as a temporary value that
%acocunts for meling in the current time step relative to the SWE output
%from the previous timestep.
%>> SNEQV  = MAX(0.,TEMP1-XM(1))

%overall this model closely replicates the reference simulation:

% % % %define default param settings:
% % % BATS_VIS_NEW = 0.95;
% % % BATS_NIR_NEW = 0.65;
% % % BATS_COSZ = 2;
% % % BATS_VIS_AGE = 0.2;
% % % BATS_NIR_AGE = 0.5;
% % % BATS_VIS_DIR = 0.4;
% % % BATS_NIR_DIR = 0.4;
% % % TAU0 = 1e6;
% % % GRAIN_GROWTH = 5000;
% % % EXTRA_GROWTH = 10;
% % % DIRT_SOOT = 0.3;
% % % SWEMX = 1;

%% load in hourly reference data:
colnames ={'DATENUM','ALBEDO','COSZ','ECAN','EDIR','ETRAN','EVB','EVC','EVG','FIRA','FSA','FSNO','FVEG','GHB','GHV','GRDFLX','HFX','IRB','IRC','IRG','LH','LWFORC','QRAIN','QSNOW','RAINRATE','SAG','SAI','SFCRNOFF','SHB','SHC','SNEQV','SNOWH','SWFORC','T2MB','T2MV','TAH','TG','TGB','TGV','TR','TRAD','TV','UGDRNOFF','SNOW_T1','SNOW_T2','SNOW_T3','ALBSND_VIZ','ALBSND_NIR','ALBSNI_VIZ','ALBSNI_NIR'};
reference_sim_data = load('/Volumes/Pruina_External_Elements/Gochis_Albedo/Data/Outputs/ReferenceHourly/reference_hourly_outputs_ObsSWDOWN.mat');
reference_sim_data = reference_sim_data.Export_Data;
ref_dates = double(reference_sim_data.DATENUM);
%% define time range to run:
start_date = datenum([2015 10 1]);
end_date = datenum([2020 9 30]);
datelist = start_date:1/24:end_date;
ndates= length(datelist);
DATEVECS = datevec(datelist);

%% load time series inputs:
%for snow age
TG = reference_sim_data.TG;
TRI_timeseries = reference_sim_data.TRI;
TRD_timeseries = reference_sim_data.TRD;
SWE = reference_sim_data.SNEQV;

%for BATS
COSZ = reference_sim_data.COSZ;

%% define range for each BATS parameter determined to be sensitive:
BATS_VIS_NEW_range = [0.8,0.95];
BATS_NIR_NEW_range = [0.55,0.80];
BATS_VIS_AGE_range = [0.6,0.99];
BATS_NIR_AGE_range = [0.6,0.99];
BATS_VIS_DIR_range = [0.2,0.8];
BATS_NIR_DIR_range = [0.2,0.8];

%% define range for each snow age parameter:
DIRTSOOT_range = [0.1 0.6];
EXTRAGROWTH_range = [5 20];
GRAINGROWTH_range = [2500 10000];
TAU0_range = [1e5 1e7];%larger TAU0 results in slower ablation degradation
SWEMX_range = [0.5 2];

rng(1)
nruns=15000;
min_ranges_p=[BATS_VIS_AGE_range(1),BATS_NIR_AGE_range(1),TAU0_range(1),BATS_NIR_DIR_range(1),BATS_VIS_NEW_range(1),BATS_NIR_NEW_range(1),DIRTSOOT_range(1),GRAINGROWTH_range(1),SWEMX_range(1)];
max_ranges_p=[BATS_VIS_AGE_range(2),BATS_NIR_AGE_range(2),TAU0_range(2),BATS_NIR_DIR_range(2),BATS_VIS_NEW_range(2),BATS_NIR_NEW_range(2),DIRTSOOT_range(2),GRAINGROWTH_range(2),SWEMX_range(2)];

[X_scaled,X_normalized]=lhsdesign_modified(nruns,min_ranges_p,max_ranges_p);

%run simulation:
ndates = length(datelist);

Store_VIZ = nan(ndates,nruns);
Store_NIR = nan(ndates,nruns);

iter=0;
for s=1:nruns
    s
    %% define input parameters:
    %for snow age
    TAU0 = X_scaled(s,3);
    GRAIN_GROWTH = X_scaled(s,8);
    EXTRA_GROWTH = 10;
    DIRT_SOOT = X_scaled(s,7);
    SWEMX = X_scaled(s,9);
    
    %for BATS
    BATS_VIS_NEW = X_scaled(s,5);
    BATS_NIR_NEW = X_scaled(s,6);
    BATS_COSZ = 2;
    BATS_VIS_AGE = X_scaled(s,1);
    BATS_NIR_AGE = X_scaled(s,2);
    BATS_VIS_DIR = 0.4;
    BATS_NIR_DIR = X_scaled(s,4);
    
    %compute snow albedo:
    [visible_snow_albedo, NIR_snow_albedo] = BATS_Snow_Age_and_Albedo(TAU0,GRAIN_GROWTH,EXTRA_GROWTH,DIRT_SOOT,SWEMX,BATS_VIS_NEW,BATS_NIR_NEW,BATS_COSZ,BATS_VIS_AGE,BATS_NIR_AGE,BATS_VIS_DIR,BATS_NIR_DIR,TG,SWE,COSZ,TRI_timeseries,TRD_timeseries);
    Store_VIZ(:,s) = visible_snow_albedo;
    Store_NIR(:,s) = NIR_snow_albedo;
    
end

%% Export
Export_Data.VIZ_albedo = Store_VIZ;
Export_Data.NIR_albedo = Store_NIR;

%export:
outdir = '/Volumes/Pruina_External_Elements/Gochis_Albedo/NoahMP_Calibration/Outputs/Small_Model_Cal/';
if exist(outdir,'dir') ==0
    CMD = ['mkdir ',outdir];
    system(CMD);
end
outfilename = 'Small_AlbMod_SnowAlb_calibration_UsingFunc.mat';
save([outdir,outfilename],'Export_Data', '-v7.3');

%% does this reproduce original code?:
original = load('/Volumes/Pruina_External_Elements/Gochis_Albedo/NoahMP_Calibration/Outputs/Small_Model_Cal/Small_AlbMod_SnowAlb_calibration.mat');
original = original.Export_Data;
original_VIZ_albedo =original.VIZ_albedo;
original_NIR_albedo =original.NIR_albedo;

IDX_VIZ = original_VIZ_albedo == Store_VIZ;
IDX_NIR = original_NIR_albedo == Store_NIR;

IDX_VIZ=IDX_VIZ(:);
IDX_NIR=IDX_NIR(:);