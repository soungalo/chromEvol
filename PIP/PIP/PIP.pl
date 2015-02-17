#!/usr/bin/perl
#-------------------------------------------------------------------------#
=begin

Program name: 	Ploidy Inference Pipeline
Version:		1.1.3
Last change:	Feb 2015

Author: 		Lior Glick
E-mail: 		liorglic@mail.tau.ac.il

Description:	This is the source code for the ploidy inference pipeline.
				The script uses the ChromEvol program in order to infer 
				extant taxa as diploid or polyploid relative to the base 
				chromosome number of the group examined.
				
Usage:			For usage instructions and examples, see the manual:
				http://www.tau.ac.il/~itaymay/cp/chromEvol/chromEvol_v2.0_manual.pdf
				
Requirements:	The script is intended for Unix\Linux machines. Make sure
				you have a working chromEvol executable, as well as the
				chromevol.pm package and the Perl packages detailed in the manual.

=cut
#-------------------------------------------------------------------------#

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin";
use chromevol;
use File::Basename;
use FindBin;
use lib "$FindBin::Bin";


my $usage = "
Usage: Ploidy Inference Pipeline
	Use control file to specify parameters
\n";
die $usage unless (scalar(@ARGV) == 1);



# Global parameters
my ($controlFile) = @ARGV;
my @paramsList = @{chromevol::read_input($controlFile)};
my ($inConsenTree,$countsFile,$inTreesFile, $infTreeNum, $simNum,
	 $paramTemplateDir, $runModels, $excludeModels, $baseNum, $ploidyCallType,
	 $duplModelsOnly, $workdir, $chromEvolExe, $plot, $cpus, $logFile) = @paramsList;


### Check input
if(check_input($countsFile,$inConsenTree,$inTreesFile)==1){
	die "ERROR: Bad input";
}
elsif(check_input($countsFile,$inConsenTree,$inTreesFile)==1){
	die "ERROR: Bad input";
}

$cpus = detect_cpus() if $cpus == 0;
if (!-d $workdir){
	system("mkdir -p $workdir") == 0
		or die "ERROR: Could not create output directory $workdir";
}
chromevol::print_params_to_log(\@paramsList,$logFile);
open (LOG, ">> $logFile") or die "ERROR: Can't write to log file '$logFile'\n";

### Run proccess
main();

### Subroutines

sub main{
	### 1. Edit counts file
	print LOG "Editing chromosome counts file...\n";
	my $countsEdit = "$workdir/".basename($countsFile)."_edit";
	chromevol::edit_counts_file($inConsenTree, $countsFile, $countsEdit) or die "ERROR: Failed to edit chromosome counts file\n$!";
	if (chromevol::check_counts($countsEdit)){
		print LOG "Error in edited counts file!\n";
		die "Error in edited counts file!\n";
	}
	
	### 2. Prepare initial chromEvol run
	print LOG "Preparing for initial chromEvol run on guide-tree...\n";
	chromevol::prepare_initial_run($inConsenTree, $countsEdit, $paramTemplateDir, $runModels, $excludeModels, $baseNum, $workdir)
		or die "ERROR: Failed when trying to prepare initial run for $inConsenTree";
	
	### 3. Run ChromEvol
	print LOG "Running chromEvol to determine best chromosome number evolution model...\n";
	opendir( D, $workdir ) or die "can't open folder $workdir for chromEvol run";
	my @chromEvolRunFiles = grep { /param_.+\.txt/ } readdir(D);
	foreach my $runFile (@chromEvolRunFiles){
		$runFile = "$workdir/$runFile";
	}
	chromevol::run_chromevol([@chromEvolRunFiles], $cpus, $chromEvolExe) or die "ERROR: Failed to complete initial chromEvol run";
	die "Initial run was not completed successfully!\n" unless chromevol::check_runs(\@chromEvolRunFiles,$workdir);
	
	### 4. Parse ChromEvol results
	print LOG "Parsing initial chromEvol run results...\n";
	# calculate max chrom number = (max in data + 10)*2
	my @chromNums = @{chromevol::calc_min_max_chrom_num($countsEdit)};
	my $maxChromNum = $chromNums[2];	# (max in data + 10)*2
	# parse results
	my $CEresultSum = "$workdir/result_sum.txt";
	my $rootFreqFile = "$workdir/root_freq.txt";
	chromevol::parse_chromevol_results($workdir, $CEresultSum, $maxChromNum, $rootFreqFile)
		or die "ERROR: Failed to parse ChromEvol results in directory $workdir";
	
	# terminate if no dupl model was chosen
	my $bestModel = chromevol::get_best_model($CEresultSum);
	if (($bestModel eq "CONST_RATE_NO_DUPL" or $bestModel eq "LINEAR_RATE_NO_DUPL") and $duplModelsOnly == 1){
		print LOG "Model $bestModel was chosen. Therefore, no polyploids can be found. Analysis terminated.\n";
		close LOG;
		open (PL, "> $workdir/ploidy.txt");
		print PL "# Ploidy inference for model $bestModel\n# No inference performed for NO_DUPL models\n";
		close PL;
		die "Model $bestModel was chosen. Therefore, no polyploids can be found. Analysis terminated.\n";
	}
		
	
	if ($inTreesFile ne "no"){
		### 5. Prepare for ChromEvol inference from real trees
		print LOG "Preparing for chromEvol inference from $infTreeNum tree topologies...\n";
		## choose trees
		my $treesListRef = chromevol::choose_trees($inTreesFile,$infTreeNum);
		if (not $treesListRef){	# not enough trees in input file
			
		}
		my $chosenTreesFile = "$workdir/chosen_trees";
		open (my $ofh, "> $chosenTreesFile") or die "ERROR: Can't write chosen trees to file\n";
		print $ofh join("\n",@{$treesListRef});
		close $ofh;
		## prepare inference
		my $inferParamTemplate = "$paramTemplateDir/param_INFER";
		my $realInferDir = "$workdir/infer";
		system("mkdir -p $realInferDir") == 0
			or die "ERROR: Failed to create inference directory $realInferDir";	
		my $realInferRunFilesRef = chromevol::prepare_infer_from_real($inTreesFile,$countsEdit, $CEresultSum,
			 $inferParamTemplate, $treesListRef, $realInferDir, $chromEvolExe)
			or die "ERROR: Failed to prepare inference from trees in file $inTreesFile";
		
		### 6. Run ChromEvol inference from real trees
		print LOG "Inferring ploidy from $infTreeNum tree topologies...\n";
		chromevol::run_chromevol($realInferRunFilesRef, $cpus, $chromEvolExe)
		 	or die "ERROR: Failed to complete chromEvol inference from real trees";
		
	}
	
	### RUN SIMULATIONS
	
	my $simulationDir = "$workdir/simulation"; # parameter and tree files
	system("mkdir -p $simulationDir") == 0
		or die "ERROR: Failed to create simulation directory $simulationDir";	

	### 7. Prepare simulations
	print LOG "Preparing to create $simNum simulation ...\n";
	my $simRunFilesListRef;
	my $simParamTemplate = "$paramTemplateDir/param_SIM";
	my $minInData = $chromNums[0];
	my $maxInData = $chromNums[1];
	my $rangeInData = $maxInData - $minInData;
	my $realInferDir = "$workdir/infer";
	if ($inTreesFile ne "no"){	# if different tree topologies are available
		$simRunFilesListRef = chromevol::prepare_simulations($realInferDir, $rootFreqFile, $simParamTemplate, $simNum,
		$maxChromNum, $rangeInData, $simulationDir, $chromEvolExe,0,0)
			or die "Error: Failed to prepare simulations\n";
	}
	else{	# if only con/MAP tree is available
		$simRunFilesListRef = chromevol::prepare_simulations($realInferDir, $rootFreqFile, $simParamTemplate, $simNum,
		$maxChromNum, $rangeInData, $simulationDir, $chromEvolExe,0,$inConsenTree)
			or die "Error: Failed to prepare simulations\n";		
	}
	
	### 8. Run simulations
	print LOG "Simulating ...\n";
	chromevol::run_chromevol($simRunFilesListRef, $cpus, $chromEvolExe)
		 or die "ERROR: Failed to complete chromEvol simulations";
	
	print LOG "Testing simulations...\n";
	my $goodSimsListRef = chromevol::test_simulations($simNum, $countsFile, $simulationDir);
	while(!$goodSimsListRef){	# while not enough good simulations have been acquired
		print LOG "$simNum reliable simulations could not be acquired. Preparing alternative simulations...\n";
		my $fakeRootFreqFile = "$workdir/fake_root_freq.txt";
		my $fakeHighestFreq = chromevol::prepare_fake_root_freq($rootFreqFile, $fakeRootFreqFile);
		print LOG "Setting fake root freq to $fakeHighestFreq\n";
		if ($inTreesFile ne "no"){
			$simRunFilesListRef = chromevol::prepare_simulations($realInferDir, $fakeRootFreqFile, $simParamTemplate, $simNum,
			$maxChromNum, $rangeInData, $simulationDir, $chromEvolExe,1,0);
		}
		else{
			$simRunFilesListRef = chromevol::prepare_simulations($realInferDir, $fakeRootFreqFile, $simParamTemplate, $simNum,
			$maxChromNum, $rangeInData, $simulationDir, $chromEvolExe,1,$inConsenTree);
			
		}		
		chromevol::run_chromevol($simRunFilesListRef, $cpus, $chromEvolExe)
		 	or die "ERROR: Failed to complete chromEvol inference from real trees";
		$rootFreqFile = $fakeRootFreqFile;
		$goodSimsListRef = chromevol::test_simulations($simNum, $countsFile, $simulationDir);
	}
	chromevol::print_good_sims($goodSimsListRef,$workdir);
	
	### 9. Parse simulations results
	my $simDataDir = "$simulationDir/simdata";
	system("mkdir -p $simDataDir") == 0
		or die "ERROR: Could not create simulation data directory $simDataDir";
	chromevol::parse_sim_results($countsEdit, $simulationDir, $goodSimsListRef, $simDataDir)
		or die "ERROR: Failed to parse simulations results in directory $simulationDir";
	
	### 10. Prepare for ChromEvol inference from simulations
	print LOG "Preparing for inference from simulated data...\n";
	my $simInferDir = "$simulationDir/infer";
	system("mkdir -p $simInferDir");
	my $inferParamTemplate = "$paramTemplateDir/param_INFER";
	my $simInferRunFilesRef = chromevol::prepare_infer_from_sim($simulationDir, $simDataDir,
	 $CEresultSum, $inferParamTemplate, $goodSimsListRef, $baseNum, $countsEdit, $simInferDir, $chromEvolExe)
	or die "ERROR: Failed to prepare inference from simulated data";
	 
	### 11. Run ChromEvol inference from simulations
	print LOG "Inferring ploidy from simulated data...\n";
	chromevol::run_chromevol($simInferRunFilesRef, $cpus, $chromEvolExe)
		 or die "ERROR: Failed to complete chromEvol inference from simulated data";
	
	### 12. Determine threshold for calling diploids/polyploids (using simulations)
	print LOG "Assessing reliability of inference...\n";
	my $thresFilePP = "$workdir/thresholds_PP";
	my $thresFileDP = "$workdir/thresholds_DP";
	my ($thresPP, $thresDP) = @{chromevol::determine_thresholds($goodSimsListRef, $simInferDir, $simulationDir,
		$thresFilePP, $thresFileDP, $ploidyCallType)}
		or die "ERROR: Failed to determine thresholds ";
	
	### 13. Compute reliability from simulations inference
	my $relSim = "$workdir/reliability_sims.txt";
	chromevol::compute_reliability_from_sim($bestModel, $workdir, $goodSimsListRef, $simInferDir,  
		$simulationDir, $thresPP, $thresDP, $relSim, $ploidyCallType)
		or die "ERROR: Failed to compute reliability from simulations in directory $simulationDir";
		
	### 14. Compute reliability from real inference
	my $relInfer = "$workdir/reliability_infer.txt";
	if ($inTreesFile ne "no"){
		my $realInferDir = "$workdir/infer";
		chromevol::compute_reliability_from_real($realInferDir, $infTreeNum, $thresPP, $thresDP, $relInfer, $ploidyCallType)
			or die "ERROR: Failed to compute reliability from real inferences in directory $realInferDir";
	}
	else{
		my $bestModelExpFile = "$workdir/$bestModel/expectations.txt";
		chromevol::compute_reliability_from_single($bestModelExpFile, $thresPP, $thresDP, $relInfer, $ploidyCallType)
			or die "ERROR: Failed to compute reliability from single inference $bestModelExpFile";
	}
	
	### 15. Summarize reliability
	my $finalOut = "$workdir/ploidy.txt";
	chromevol::summarize_reliability($relSim,$relInfer,$bestModel,$countsEdit,$finalOut)
		or die "ERROR: Failed to summarize reliability";
	
	### 16. create ancestral trees files and parameter distributions
	chromevol::get_ancestral_trees($workdir);
	if ($inTreesFile ne "no"){
		my $realInferDir = "$workdir/infer";	
		my $paramDistFile = "$workdir/param_distribution";
		chromevol::get_params_distribution($realInferDir,$infTreeNum,$paramDistFile)
			or die "ERROR: Failed to summarize parameter distributions!\n";
	}
	
	print LOG "# Process completed successfuly";
}
close LOG;

sub check_input{
	my ($countsFile, $conTreeFile, $treesFile) = @_;
	my $inputError = 0;
	if (chromevol::check_counts($countsFile)==1){
		$inputError=1;
	}
	if (chromevol::check_con_tree($conTreeFile)==1){
		$inputError=1;
	}
	if ($treesFile ne "no" && chromevol::check_tree_sample($treesFile)==1){
		$inputError=1;
	}
	return $inputError;
}

sub detect_cpus{
	my $cpusString = `grep processor /proc/cpuinfo`;
	my @procLines = split(/\n/,$cpusString);
	my $cpus = scalar(@procLines);
	return $cpus;
}


