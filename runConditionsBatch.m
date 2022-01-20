clear

%%% For 2 micron data
% % Varying RI
dataFolder = '/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/AnalysisVolumes/2 micron/Varying RI Profile/';
interConeValueToUse = 1.40;

% metaFile = 'Cone_1000_nm_Cone_0_SD_GRIN_cylinder.mat'; graded_RI_trace_and_test
% metaFile = 'Cone_1000_nm_Cone_0_SD_GRIN_cylinder_m36.mat'; graded_RI_trace_and_test
% metaFile = 'Cone_1000_nm_Cone_0_SD_GRIN_cylinder_p36.mat'; graded_RI_trace_and_test
% metaFile = 'Cone_1000_nm_Cone_0_SD_GRIN_cylinder_normalised_max.mat'; graded_RI_trace_and_test
metaFile = 'Cone_1000_nm_Cone_0_SD_GRIN_cylinder_normalised_80.mat'; graded_RI_trace_and_test
% metaFile = 'Cylinder_1000_nm_Cone_0_SD_GRIN_cylinder.mat'; graded_RI_trace_and_test

% metaFile = 'Cone_1000_nm_Cone_0_SD_GRIN_both.mat'; graded_RI_trace_and_test
% metaFile = 'Cone_1000_nm_Cone_0_SD_GRIN_linear.mat'; graded_RI_trace_and_test
% metaFile = 'Cone_1000_nm_Cone_0_SD_GRIN_radialTop.mat'; graded_RI_trace_and_test
% metaFile = 'Cone_1000_nm_Cone_0_SD_GRIN_radialBase.mat'; graded_RI_trace_and_test
% metaFile = 'Cone_1000_nm_Cone_0_SD_GRIN_radialTop_TipCorrection.mat'; graded_RI_trace_and_test

% for interConeValueToUse = 1.42:0.02:1.48
%     metaFile = 'Cone_1000_nm_Cone_0_SD_GRIN_radialTop.mat'; graded_RI_trace_and_test
% 
%     metaFile = 'Cone_1000_nm_Cone_0_SD_GRIN_cylinder.mat'; graded_RI_trace_and_test
% end

% internal structures
% dataFolder = '/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/AnalysisVolumes/2 micron/Extra internal structures/Radial RI/';
% interConeValueToUse = 1.40;

% metaFile = 'Cone_CinC_1000_nm_Cone_0_SD_GRIN_radialTop.mat'; graded_RI_trace_and_test
% metaFile = 'Cone_CinC_EC_1000_nm_Cone_0_SD_GRIN_radialTop.mat'; graded_RI_trace_and_test
% metaFile = 'Cone_EC_1000_nm_Cone_0_SD_GRIN_radialTop.mat'; graded_RI_trace_and_test

% % Uniform
% % Varying shape SD
% dataFolder = '/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/AnalysisVolumes/2 micron/Varying Cone SD/';
% interConeValueToUse = 1.40;
 
% metaFile = 'Cone_1000_nm_Cone_-4_SD_Uniform_1.52.mat'; graded_RI_trace_and_test
% metaFile = 'Cone_1000_nm_Cone_-2_SD_Uniform_1.52.mat'; graded_RI_trace_and_test
% metaFile = 'Cone_1000_nm_Cone_0_SD_Uniform_1.52.mat'; graded_RI_trace_and_test
% metaFile = 'Cone_1000_nm_Cone_2_SD_Uniform_1.52.mat'; graded_RI_trace_and_test
% metaFile = 'Cone_1000_nm_Cone_4_SD_Uniform_1.52.mat'; graded_RI_trace_and_test
 
% for interConeValueToUse = 1.42:0.02:1.48
%     metaFile = 'Cone_1000_nm_Cone_0_SD_Uniform_1.52.mat'; graded_RI_trace_and_test
% end

% % internal structures
% dataFolder = '/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/AnalysisVolumes/2 micron/Extra internal structures/Uniform RI/';
% interConeValueToUse = 1.40;

% metaFile = 'Cone_CinC_1000_nm_Cone_0_SD_Uniform_1.52.mat'; graded_RI_trace_and_test
% metaFile = 'Cone_CinC_EC_1000_nm_Cone_0_SD_Uniform_1.52.mat'; graded_RI_trace_and_test
% metaFile = 'Cone_EC_1000_nm_Cone_0_SD_Uniform_1.52.mat'; graded_RI_trace_and_test

close all

