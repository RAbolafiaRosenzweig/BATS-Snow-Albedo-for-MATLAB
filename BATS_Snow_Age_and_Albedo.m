function [visible_snow_albedo, NIR_snow_albedo] = BATS_Snow_Age_and_Albedo(TAU0,GRAIN_GROWTH,EXTRA_GROWTH,DIRT_SOOT,SWEMX,BATS_VIS_NEW,BATS_NIR_NEW,BATS_COSZ,BATS_VIS_AGE,BATS_NIR_AGE,BATS_VIS_DIR,BATS_NIR_DIR,TG,SWE,COSZ,TRI_timeseries,TRD_timeseries)

%%================================================================================
%% defintions of input parameter definitions for snow age (FAGE) calculation
%%================================================================================
%TAU0: parameter from from Yang97 eqn. 10a, default = 1.e6
%GRAIN_GROWTH: growth from vapor diffusion Yang97 eqn. 10b, default = 5000
%EXTRA_GROWTH: extra growth near freezing Yang97 eqn. 10c, default = 10
%DIRT_SOOT: dirt and soot term Yang97 eqn. 10d, default = 0.3
%SWEMX: new snow mass to fully cover old snow (mm), default = 1.00

%%================================================================================
%% defintions of input parameter definitions for snow albedo calculation
%%================================================================================
%BATS_VIS_NEW: new snow visible albedo
%BATS_NIR_NEW: new snow NIR albedo
%BATS_COSZ: zenith angle snow albedo adjustment; b in Yang97 eqn. 15
%BATS_VIS_AGE: age factor for diffuse visible snow albedo Yang97 eqn. 17
%BATS_NIR_AGE: age factor for diffuse NIR snow albedo Yang97 eqn. 18
%BATS_VIS_DIR: cosz factor for direct visible snow albedo Yang97 eqn. 15
%BATS_NIR_DIR: cosz factor for direct NIR snow albedo Yang97 eqn. 16

%%================================================================================
%% defintions of input time series required to compute snow albedo
%%================================================================================
%TG: ground temperature (k) (from Noah-MP reference simulation)
%SWE: snow water equivalent (mm) (from Noah-MP reference simulation)
%COSZ: cosine solar zenith angle (0-1) (from Noah-MP reference simulation)
%TRI_timeseries: transmitted solar radiation: diffuse (w/m2) (from Noah-MP reference simulation)
%TRD_timeseries: transmitted solar radiation: direct (w/m2)(from Noah-MP reference simulation)

%%================================================================================
%% define other fixed parameters
%%================================================================================
%define number of bands = 2 (visible & NIR)
NBAND = 2;
%define duration of hourly timestep:
DT = 3600; %seconds
%define freezing temperature of water
TFRZ = 273.16; %as hardcoded in Noah-MP

%%================================================================================%%================================================================================
%start time loop:
iter = 0;

ndates = length(SWE);
visible_snow_albedo =nan(ndates,1);
NIR_snow_albedo =nan(ndates,1);
for timestep = 1:ndates
    %%================================================================================
    %% perform BATS snow age calculation
    %%================================================================================
    if COSZ(timestep) > 0 && SWE(timestep) > 0
        iter = iter+1;
        if iter==1
            TAUSS = 0;
        end
        DELA0 = DT/TAU0;
        ARG   = GRAIN_GROWTH*(1./TFRZ-1./TG(timestep));
        AGE1  = exp(ARG);
        AGE2  = exp(min(0,EXTRA_GROWTH*ARG));
        AGE3  = DIRT_SOOT;
        TAGE  = AGE1+AGE2+AGE3;
        DELA  = DELA0*TAGE;
        SNEQV=SWE(timestep);
        if timestep == 1
            SNEQVO = 0;
        else
            SNEQVO=SWE(timestep-1);
            if SNEQVO <=1e-2
                SNEQVO=0;
            end
        end
        DELS  = max(0.0,SNEQV-SNEQVO) / SWEMX;
        SGE   = (TAUSS+DELA)*(1.0-DELS);
        TAUSS = max(0,SGE);
        FAGE= TAUSS/(TAUSS+1);
        
        %%================================================================================
        %% perform BATS snow albedo calculations
        %%================================================================================
        %initialize
        ALBSND(1: NBAND) = 0;
        ALBSNI(1: NBAND) = 0;
        
        %compute:
        SL=BATS_COSZ;
        SL1=1./SL;
        SL2=2.*SL;
        CF1= ((1+SL1)/(1+SL2*COSZ(timestep))-SL1);
        FZEN=max(CF1,0);
        
        %diffuse
        ALBSNI(1)=BATS_VIS_NEW*(1-BATS_VIS_AGE*FAGE);
        ALBSNI(2)=BATS_NIR_NEW*(1-BATS_NIR_AGE*FAGE);
        %direct
        ALBSND(1)=ALBSNI(1)+BATS_VIS_DIR*FZEN*(1-ALBSNI(1));
        ALBSND(2)=ALBSNI(2)+BATS_NIR_DIR*FZEN*(1-ALBSNI(2));
        
        TRI = TRI_timeseries(timestep);
        TRD = TRD_timeseries(timestep);
        
        if TRI==0 && TRD==0
            Weight_I = 0;
            Weight_D = 0;
        else
            Weight_I = TRI/(TRI+TRD);
            Weight_D = TRD/(TRI+TRD);
        end
        
        visible_snow_albedo(timestep) = ALBSND(1)*Weight_D + ALBSNI(1)*Weight_I;
        NIR_snow_albedo(timestep) = ALBSND(2)*Weight_D + ALBSNI(2)*Weight_I;
        
        assert(visible_snow_albedo(timestep)>=0 & visible_snow_albedo(timestep)<=1,'albedo is out of range');
        assert(NIR_snow_albedo(timestep)>=0 & NIR_snow_albedo(timestep)<=1,'albedo is out of range');
    else %cosz<0, SWDOWN = 0
        visible_snow_albedo(timestep) =0;
        NIR_snow_albedo(timestep)=0;
    end
end