# nbTransmission-updates
Folder nbTransmission contains updated versions of performNB, runCV, and nbProbabilities (from package nbTransmission) to generate adjusted odds ratios via a bootstrapped logisitic regression model. nbProbabilities_BS performs a full bootstrap and nbProbabilities_sBS performs a modified bootstrap that uses all cases and a smaller sample of controls. 
<br />
<br />
Folder nbSimulation contains code to simulate outbreaks and covariates called in SimRunEst.R. SimCovariatesNBEffect_corr.R simulates covariates with two covariates correlated in order to assess the adjusted odds ratios. 
<br />
<br />
PerformSimulationEst_corr calls SimRunEst to run simulations of the original methods in R package nbTransmission as well as bootstrapped and non-bootstrapped logisitic regression estimated adjusted odds ratios within the iterative framework. 
