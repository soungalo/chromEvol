package chromevol;

use strict;
use warnings;
use Bio::TreeIO;
use List::Util qw(shuffle);
use Parallel::ForkManager;
use POSIX;
use List::Util qw( min max );
use File::Copy;
use File::Basename;

# sub edit counts file
# ensures the propper format of the counts file as input for chromEvol
sub edit_counts_file{
	my ($consenTree, $countsFile, $outCounts) = @_;
	
	my @taxaInTree = @{get_tree_taxa($consenTree)};
	my %taxaInCounts = %{parse_counts_file($countsFile)};
	
	my %countsHash;
	foreach my $treeTaxon (@taxaInTree){
		if (exists $taxaInCounts{$treeTaxon}){
			$countsHash{$treeTaxon} = $taxaInCounts{$treeTaxon};
		}
		else{
			$countsHash{$treeTaxon} = "x";
		}
	}
	
	open (OUT, "> $outCounts") or die "Can't write edited counts file $outCounts";
	my @keys = keys(%countsHash);
	foreach my $key (@keys){
		print OUT ">$key\n".$countsHash{$key}."\n";
	}
	close(OUT);
	return 1;
}

# sub get trees taxa
# get an array of the taxa names in the input newick tree
sub get_tree_taxa{
	my ($treesFile) = @_;
	my ($format) = qw(newick);
	my $in = Bio::TreeIO->new(-file => $treesFile, -format => $format);

	my $t = $in->next_tree;
	my @treeTaxa;
	my @taxaObject = $t->get_leaf_nodes;
	for (my $i = 0; $i<=$#taxaObject; $i++) {
		my $taxon = $taxaObject[$i]{"_id"};
		push(@treeTaxa,$taxon) unless not defined $taxon;
	}
	return [@treeTaxa];
}

# sub parse counts file
# returns a hash with the counts file names and counts
sub parse_counts_file{
	my $fastaFile = $_[0];
	my %counts;
	open (IN, "< $fastaFile") or die "Can't open $fastaFile";
	my $fastaLine = <IN>;
	while (defined $fastaLine){
		chomp $fastaLine;
		if ($fastaLine =~ m/^>(.+)/){
			my $taxon = $1;
			$taxon =~ s/\s/_/;	# in case taxon name contains space, replace with '_'
			$fastaLine = <IN>;
			chomp $fastaLine;
			$counts{$taxon} = $fastaLine;
		}
		$fastaLine = <IN>;
	}
	close IN;
	return {%counts};
}

# sub check counts
# checks that the input counts file is in the expected format, and raises a warning otherwise
sub check_counts{
	my $countsFile = $_[0];
	my $error = 0;	# 0 if file is OK, 1 if error
	open (CNT , "< $countsFile") or die "ERROR: Can't open chromosome counts file $countsFile\n$!";
	my $line = <CNT>;
	while (defined $line){
		chomp $line;
		if ($line =~ m/^>([^\s\t]+)/){	# check that taxon name starts with > and does not include spaces
			my $taxon = $1;			
			$line = <CNT>;
			chomp $line;
			# check that count is either a number or an 'x'
			if ($line =~ m/[1-9]\d*/ or $line eq "x" or $line eq "X"){
				
			}
			else{
				warn "Bad chromosome count in line '$line' for taxon $taxon\n";
				$error = 1;
				last;
			}
		}
		else{
			warn "Bad taxon name in line '$line'\n";
			$error = 1;
			last;
		}
		$line = <CNT>;
	}
	# check that file has at least one count
	my %countsHash = %{parse_counts_file($countsFile)};
	my @counts = values(%countsHash);
	@counts = grep(($_ ne "x" and $_ ne "X"), @counts);
	unless (scalar(@counts) >= 1){
		$error = 1;
		warn "No counts found in counts file\n";
	}
	
	return $error;
}

# sub check consensus tree
# checks that the input consensus tree is in the expected format, and raises a warning otherwise
sub check_con_tree{
	my $conTreeFile = $_[0];
	my $error = 0;	# 0 if file is OK, 1 if error
	open (CON, "< $conTreeFile") or die "ERROR: Can't open consensus tree file $conTreeFile\n$!";
	my @lines = <CON>;
	chomp @lines;
	if (scalar(@lines)>1){
		warn "More than one tree in consensus tree file!\n";
		$error = 1;
	}
	elsif (scalar(@lines)==0){
		warn "Consensus tree file is empty!\n";
		$error = 1;
	}
	else{
		# CHECK IF TREE HAS BOOTSTRAP
		my ($format) = qw(newick);
		my $in = Bio::TreeIO->new(-file => $conTreeFile, -format => $format);
		my $tree = $in->next_tree;
		my @internalNodes =  grep { ! $_->is_Leaf } $tree->get_nodes;
		foreach my $node (@internalNodes){
			my $bs = $node->bootstrap;
			my $id = $node->id;
			if ($bs or $id){
				$error = 1;
				warn "Consensus tree has bootstrap values!\n";
				last;
			}
		}	
	}
	return $error;
}

# sub check trees sample
# checks that the input trees sample file is in the expected format, and raises a warning otherwise
sub check_tree_sample{
	my $treesFile = $_[0];
	my $error = 0;	# 0 if file is OK, 1 if error
	open (CON, "< $treesFile") or die "ERROR: Can't open trees sample file $treesFile\n$!";
	my @lines = <CON>;
	chomp @lines;
	if (scalar(@lines) == 1){
		warn "Only one tree in trees file!\n";
		$error = 1;
	}
	elsif (scalar(@lines)==0){
		warn "Trees file is empty!\n";
		$error = 1;
	}
	else{
		# CHECK IF TREE HAS BOOTSTRAP
		my ($format) = qw(newick);
		my $in = Bio::TreeIO->new(-file => $treesFile, -format => $format);
		my $tree = $in->next_tree;
		my @internalNodes =  grep { ! $_->is_Leaf } $tree->get_nodes;
		foreach my $node (@internalNodes){
			my $bs = $node->bootstrap;
			my $id = $node->id;
			if ($bs or $id){
				$error = 1;
				warn "One or more trees has bootstrap values!\n";
				last;
			}
		}	
	}
	return $error;
}

# sub prepare initial chromevol run
# prepares the required files for the initial chromevol run
sub prepare_initial_run{
	### file and dirs
	my ($treeFile, $countFile, $paramDir, $runModels, $excludeModels, $baseNum, $outDir) = @_;
	my @modelsList;
	if ($runModels eq "ALL"){
		@modelsList = ("CONST_RATE","CONST_RATE_DEMI","CONST_RATE_DEMI_EST","CONST_RATE_NO_DUPL",
		"LINEAR_RATE","LINEAR_RATE_DEMI","LINEAR_RATE_DEMI_EST","LINEAR_RATE_NO_DUPL",
		"BASE_NUM","BASE_NUM_DUPL");
	}
	elsif ($runModels eq "CONST"){
		@modelsList = ("CONST_RATE","CONST_RATE_DEMI","CONST_RATE_DEMI_EST","CONST_RATE_NO_DUPL",
		"BASE_NUM","BASE_NUM_DUPL");		
	}
	elsif ($runModels eq "LINEAR"){
		@modelsList = ("LINEAR_RATE","LINEAR_RATE_DEMI","LINEAR_RATE_DEMI_EST","LINEAR_RATE_NO_DUPL");		
	}
	elsif ($runModels eq "NO_BASE"){
		@modelsList = ("CONST_RATE","CONST_RATE_DEMI","CONST_RATE_DEMI_EST","CONST_RATE_NO_DUPL",
		"LINEAR_RATE","LINEAR_RATE_DEMI","LINEAR_RATE_DEMI_EST","LINEAR_RATE_NO_DUPL");
	}
	else{
		@modelsList = split(" ",$runModels);
	}
	if ($excludeModels != 0){
		my @modelsToExclude = split(" ",$excludeModels);
		foreach my $exModel (@modelsToExclude){
			@modelsList = grep{$exModel ne $_} @modelsList;	# substract @modelsList - @modelsToExclude
		}		
	}
	
	### parameter files
	my @paramFiles;
	foreach my $model (@modelsList){
		push  (@paramFiles, 'param_'.$model);
	}
	### ChromEvol run parameters
	### NB: NOT EXHAUSTIVE LIST OF PARAMETERS
	### 	SEE <http://www.tau.ac.il/~itaymay/cp/chromEvol/> FOR DETAILS
	my $maxChrNum = -1;	# highest count in file + 1
	my $minChrNum = -1;	# lowest count in file - 1
	
	### RUN
	foreach my $pFile (@paramFiles){
		open OUT, "> $outDir/$pFile.txt" or die "ERROR: failed to create $outDir/$pFile...";
		print OUT &parseParamFile($paramDir."$pFile",$outDir,$treeFile);
		close OUT;
	}
	
	
	### SUBROUTINES
	sub parseParamFile{
		no warnings;
		open FILE, $_[0] or die "ERROR: failed to open $_[0] ...";
		my $outDir = $_[1];
		my $treeFile = $_[2];
		undef $/;
		$_ = <FILE>;
		close FILE;
		$/ = "\n";
		s/\<OUT\_DIR\>/$outDir/g;
		s/\<TREE\_FILE\>/$treeFile/;
		s/\<CNT\_FILE\>/$countFile/;
		s/\<MAX\_CHR\_NUM\>/$maxChrNum/;
		s/\<MIN\_CHR\_NUM\>/$minChrNum/;
		if ($baseNum != 0){	# if user specifies base number
			s/<BASE_NUM>/$baseNum/;
			s/<OPT_BASE_NUM>/0/;	# don't optimize base number
		}
		else{
			# find min count in data - will be used as start point for base number optimization
			my @minMaxCounts = @{calc_min_max_chrom_num($countFile,4)};
			my $minInData = $minMaxCounts[0];
			if (not $minInData){
				return 0;		# raise error if there are no counts >= 4
			}

			s/<BASE_NUM>/$minInData/;
			s/<OPT_BASE_NUM>/1/;	# optimize base number
		}
		$_ .= "\n";
		$/ = "\n";
		return $_;
	}
	use warnings;
	return 1;
}


# sub calculate max chromosome number
# calculates the maximun chromosome number allowed in later inferences = (max in data + 10)*2
sub calc_min_max_chrom_num{
	my $countsEdit = $_[0];
	my $lower_bound = $_[1];
	
	my %origCountsHash = %{parse_counts_file($countsEdit)};
	my @counts = values(%origCountsHash);
	my @numericCounts;	# will include all counts in file, but not 'x'
	foreach my $c (@counts){
		if ($c =~ m/^\d+$/){
			push(@numericCounts, $c);
		}
		elsif ($c =~ m/=/){	# intraspecific variation
			my @variations = split(/_/,$c);
			foreach my $v (@variations){
				my @parts = split(/=/,$v);
				push (@numericCounts,$parts[0]) if ($parts[0] =~ m/^\d+$/)
			}
		}
	}
	
	
	if ($lower_bound){
		@numericCounts = grep {$_ >= $lower_bound} @numericCounts;
	}
	
	my @sortCounts = sort {$b <=> $a} @numericCounts;
	my $maxInData = $sortCounts[0];
	my $minInData = $sortCounts[-1];
	my $maxChromNum = ($maxInData + 10)*2;
	my @returnArr = ($minInData,$maxInData,$maxChromNum);
	return \@returnArr;
}

# sub parse chromevol results
# parse the output of initial run to find best model
sub parse_chromevol_results{
	### FIELD
	my ($inDir, $outFile, $maxChrNbrSim, $freqFile) = @_;
	
	### RUN

		no warnings;
		my @models = ('CONST_RATE', 'CONST_RATE_DEMI', 'CONST_RATE_DEMI_EST', 'CONST_RATE_NO_DUPL', 'LINEAR_RATE', 'LINEAR_RATE_DEMI', 'LINEAR_RATE_DEMI_EST', 'LINEAR_RATE_NO_DUPL','BASE_NUM', 'BASE_NUM_DUPL');
		my @modelInfo = ();
		foreach(@models){
			my $resFile = "$inDir/$_/chromEvol.res";
			if(-e $resFile){
				push @modelInfo, &parseChromEvolResults($resFile, $_);
			}else{
				push @modelInfo, 'NA';
			}
		}
		my ($MODEL_INFO, $BEST_MODEL, $BEST_FREQ) = &getBestModel(@modelInfo);
		### print ChromEvol results to file
		open MODEL, ">$outFile" or die "ERROR: failed to create $outFile...";
		$MODEL_INFO =~ s/:/\t/g;
		print MODEL "$MODEL_INFO\n";
		close MODEL;
		### create root chromosome number frequency file
		open FREQ, ">$freqFile" or die "ERROR: failed to create $freqFile...";
		my @freqs = split m/\s+/, $BEST_FREQ;
		$freqs[0] =~ m/F\[(\d+)\]/;
		my $lowestRoot = $1;
		$freqs[(scalar(@freqs)-1)] =~ m/F\[(\d+)\]/;
		my $highestRoot = $1;
		### add to file F[R], where 1 <= R < lowestRoot
		for(my $i = 1; $i < $lowestRoot; $i++){
			print FREQ "F[$i]=0\n";
		}
		### add to file F[R], where lowestRoot <= R <= highestRoot
		for(my $i = $lowestRoot; $i <= $highestRoot; $i++){
			print FREQ "$freqs[$i-$lowestRoot]\n";
		}
		### add to file F[R], where highestRoot < R <= maxChrNbrSim
		for(my $i = $highestRoot+1; $i <= $maxChrNbrSim; $i++){
			print FREQ "F[$i]=0\n";
		}
		close FREQ;
		use warnings;
	
	
	### parse ChromEvol output results file
	sub parseChromEvolResults{
		my $file = $_[0];
		my $model = $_[1];
		my $lnLik;
		my $AIC;
		my $totTrLen;
		my $rootFreq;
		my ($gain, $gainL, $loss, $lossL, $dupl, $demi,$baseNum,$baseNumR) = ('na', 'na', 'na', 'na', 'na', 'na','na','na');
		### parse 'chromEvol.res'
		open FILE, $file or die "ERROR: failed to open $file...";
		while(<FILE>){
			chomp;
			if( m/^LogLikelihood = (.+)/ ){
				$lnLik = $1;
			}elsif( m/^AIC \(Akaike information criterion\) = (.+)/ ){
				$AIC = $1;
			}elsif( m/^GAIN_CONST\s+(.+)/ ){
				$gain = $1;
			}elsif( m/^GAIN_LINEAR\s+(.+)/ ){
				$gainL = $1;
			}elsif( m/^LOSS_CONST\s+(.+)/ ){
				$loss = $1;
			}elsif( m/^LOSS_LINEAR\s+(.+)/ ){
				$lossL = $1;
			}elsif( m/^DUPL\s+(.+)/ ){
				$dupl = $1;
			}elsif( m/^HALF_DUPL\s+(.+)/ ){
				$demi = $1;
			}elsif( m/\#total\stree\slength\s\=\s(\d+)/ ){
				$totTrLen = $1;
			}elsif( m/^F\[\d+\]/ ){
				$rootFreq = $_;
			}
			elsif (m/^BASE_NUMBER_R\s+(.+)/ ){
				$baseNumR = $1;
			}
			elsif (m/^BASE_NUMBER\s+(.+)/ ){
				$baseNum = $1;
			}
			
		}
		close FILE;
		return join(':', ($model, $lnLik, $AIC, $totTrLen, $gain, $gainL, $loss, $lossL, $dupl, $demi, $baseNum, $baseNumR, $rootFreq));
	}
	
	### identify best-fitting model based on lowest AIC score
	sub getBestModel{
		my @lines = @_;
		my @modLines = ();
		my $bestModel;
		my $bestFreq;
		my $lowestAIC;
		my @AIC = ();
		### find lowest AIC
		foreach(@lines){
			push @AIC, (split m/:/, $_)[2];
		}
		@AIC = sort {$a <=> $b} @AIC;
		$lowestAIC = $AIC[0];
		### add delta AIC
		foreach(@lines){
			next if $_ eq 'NA';
			my @tmp = split m/:/, $_;
			if($tmp[2] - $lowestAIC == 0){
				$bestModel = $tmp[0];
				$bestFreq = $tmp[12];
			}
			push @modLines, join(':', ($tmp[0], $tmp[1], $tmp[2], $tmp[2] - $lowestAIC, $tmp[3], @tmp[4 .. 11]));
		}
		return(join("\n", @modLines), $bestModel, $bestFreq);
	}
}

# sub distribute processes
# distributes a list of processes over a given number of CPUs (returns a distribution hash) 
sub distribute_proc{
	my ($procListRef, $cpus) = @_;
	my @procList = @$procListRef;	# an array of commands
	my %distribHash;
	my $cpuInd = 1;
	# create the basic data structure
	while ($cpuInd <= $cpus){
		$distribHash{"CPU$cpuInd"} = [];
		$cpuInd++;
	}
	# distribute processes to cpus
	my $procInd = 0;
	$cpuInd = 1;
	while ($procInd < scalar(@procList)){
		my @cpuProc = @{$distribHash{"CPU$cpuInd"}};
		push (@cpuProc, $procList[$procInd]);
		$distribHash{"CPU$cpuInd"} = [@cpuProc];
		$procInd++;
		$cpuInd++;
		$cpuInd = 1 if ($cpuInd > $cpus);
	}
	# remove empty CPUs
	foreach (keys(%distribHash)){
		my @cpuProcs = @{$distribHash{$_}};
		delete $distribHash{$_} if (scalar(@cpuProcs)==0);
	}
	
	return {%distribHash};
}

# sub run chromevol
# runs a list of chromevol commands in parallel CPUs
sub run_chromevol{
	my ($runFilesListRef, $numOfCpus, $chromevolExe) = @_;
	my %commandHash = %{distribute_proc($runFilesListRef, $numOfCpus)};
	my @cpus = keys(%commandHash);
	my $pm = new Parallel::ForkManager($numOfCpus);
	foreach my $cpu (@cpus){
		$pm->start and next;
		my @runFiles = @{$commandHash{$cpu}};
		foreach my $runFile(@runFiles){
			system("$chromevolExe", "$runFile") == 0
				or die "ERROR: Failed when running chromEvol on file $runFile";
		}
		$pm->finish;
	}
	$pm->wait_all_children;
	return 1;
}

# sub prepare infer from real trees
# prepares required files for inferring ploidy from different tree topologies
sub prepare_infer_from_real{
	### FIELD
	my $inTreesFile = $_[0];
	my $countsFile = $_[1];
	my $summaryFile = $_[2];
	my $paramFile = $_[3];
	my $treesListRef = $_[4];
	my $outDir = $_[5];
	my $chromEvolExe = $_[6];
	#my @trees = ();  # keep list of trees in Newick format
	
	
	### RUN
	### parse ChromEvol results summary
	# 1) model
	# 2) lnLik
	# 3) AIC
	# 4) dAIC
	open RES, $summaryFile or die "ERROR: failed to open $summaryFile...";
	my $initialParamValues = '';
	while(<RES>){
	        chomp;
	        my @tmp = split m/\t/, $_;
		### get best model
	        if($tmp[3] == 0){
			if($tmp[0] eq 'CONST_RATE'){
				$initialParamValues .= "_gainConstR $tmp[5]\n";
				$initialParamValues .= "_lossConstR $tmp[7]\n";
				$initialParamValues .= "_duplConstR $tmp[9]\n";
			}elsif($tmp[0] eq 'CONST_RATE_DEMI'){
				$initialParamValues .= "_gainConstR $tmp[5]\n";
				$initialParamValues .= "_lossConstR $tmp[7]\n";
				$initialParamValues .= "_duplConstR $tmp[9]\n";
				$initialParamValues .= " _demiPloidyR -2\n";
			}elsif($tmp[0] eq 'CONST_RATE_DEMI_EST'){
				$initialParamValues .= "_gainConstR $tmp[5]\n";
				$initialParamValues .= "_lossConstR $tmp[7]\n";
				$initialParamValues .= "_duplConstR $tmp[9]\n";
				$initialParamValues .= " _demiPloidyR $tmp[10]\n";
			}elsif($tmp[0] eq 'CONST_RATE_NO_DUPL'){
				$initialParamValues .= "_gainConstR $tmp[5]\n";
				$initialParamValues .= "_lossConstR $tmp[7]\n";
			}elsif($tmp[0] eq 'LINEAR_RATE'){
				$initialParamValues .= "_gainConstR $tmp[5]\n";
				$initialParamValues .= "_gainLinearR $tmp[6]\n";
				$initialParamValues .= "_lossConstR $tmp[7]\n";
				$initialParamValues .= "_lossLinearR $tmp[8]\n";
				$initialParamValues .= "_duplConstR $tmp[9]\n";
			}elsif($tmp[0] eq 'LINEAR_RATE_DEMI'){
				$initialParamValues .= "_gainConstR $tmp[5]\n";
				$initialParamValues .= "_gainLinearR $tmp[6]\n";
				$initialParamValues .= "_lossConstR $tmp[7]\n";
				$initialParamValues .= "_lossLinearR $tmp[8]\n";
				$initialParamValues .= "_duplConstR $tmp[9]\n";
				$initialParamValues .= "_demiPloidyR -2\n";
			}elsif($tmp[0] eq 'LINEAR_RATE_DEMI_EST'){
				$initialParamValues .= "_gainConstR $tmp[5]\n";
				$initialParamValues .= "_gainLinearR $tmp[6]\n";
				$initialParamValues .= "_lossConstR $tmp[7]\n";
				$initialParamValues .= "_lossLinearR $tmp[8]\n";
				$initialParamValues .= "_duplConstR $tmp[9]\n";
				$initialParamValues .= "_demiPloidyR $tmp[10]\n";
			}elsif($tmp[0] eq 'LINEAR_RATE_NO_DUPL'){
				$initialParamValues .= "_gainConstR $tmp[5]\n";
				$initialParamValues .= "_gainLinearR $tmp[6]\n";
				$initialParamValues .= "_lossConstR $tmp[7]\n";
				$initialParamValues .= "_lossLinearR $tmp[8]\n";
			}
			elsif($tmp[0] eq 'BASE_NUM'){
				$initialParamValues .= "_gainConstR $tmp[5]\n";
				$initialParamValues .= "_lossConstR $tmp[7]\n";
				$initialParamValues .= "_baseNumberR $tmp[12]\n";
				$initialParamValues .= "_baseNumber $tmp[11]\n";
			}
			elsif($tmp[0] eq 'BASE_NUM_DUPL'){
				$initialParamValues .= "_gainConstR $tmp[5]\n";
				$initialParamValues .= "_lossConstR $tmp[7]\n";
				$initialParamValues .= "_duplConstR $tmp[9]\n";
				$initialParamValues .= "_baseNumberR $tmp[12]\n";
				$initialParamValues .= "_baseNumber $tmp[11]\n";
			}
			
	        }
	}
	close RES;
	
	### parse template parameter file for ChromEvol
	open PARAM, $paramFile or die "ERROR: failed to open $paramFile ...";
	undef $/;
	my $paramText = <PARAM>;
	$paramText =~ s/\<PARAMETERS\>/$initialParamValues/;
	$/ = "\n";
	close PARAM;
	
	### create separate tree files
	open (TREE, "< $inTreesFile") or die "ERROR: failed to open $inTreesFile ...";
	my @treeLines = <TREE>;
	chomp @treeLines;
	@treeLines = shuffle(@treeLines);
	close TREE;
	my @treesList = @{$treesListRef};
	my $i = 1;
	foreach my $n (@treesList){
		my $tree = $treeLines[$n-1];
		my $treeFile = "$outDir/tree_$i";
		open (OUT, "> $treeFile") or die "ERROR: failed to create $treeFile...";
		print OUT "$tree";
		close OUT;
		$i++;
	}
	
	
	### print param file and create param files list
	
	my @runFilesList;
	#open CMD, ">$cmdFile" or die "ERROR: failed to open $cmdFile...";
	$i = 1;
	foreach my $n (@treesList){
        ### print param to file
        my $outInferDir = $outDir.'/dir_'.$i.'/';
        my $inTreeFile = "$outDir/tree_$i";
        my $outParamFile = $outDir.'/param_'.$i;
        my $newParamText = $paramText;
        open OUT, ">$outParamFile" or die "ERROR: failed to create $outParamFile...";
        $newParamText =~ s/\<OUT_DIR\>/$outInferDir/;
        $newParamText =~ s/\<FAC_FILE\>/$countsFile/;
        $newParamText =~ s/\<TREE_FILE\>/$inTreeFile/;
        print OUT "$newParamText\n";
        close OUT;
        ### print cmds to file
        push(@runFilesList, $outParamFile);
        $i++;
	}
	#close CMD;
	return \@runFilesList;
}


# sub check file
# checks if a file is empty
sub check_file{
	my $inFile = $_[0];
	open F, "< $inFile" or die "ERROR: Can't open file $inFile";
	my @lines = <F>;
	close F;
	
	if (scalar(@lines) == 0){
		return 0;
	}
	elsif (scalar(@lines) > 0){
		return 1;
	}
	else{
		die "ERROR: Can't open file $inFile";
	}
}

# sub test simulations
# Tests each simulation and creates a list of good simulations
sub test_simulations{
	my ($simNum,$countsFile,$allSimDir) = @_;
	
	# find chromosome range in original counts
	my @countsParse = @{calc_min_max_chrom_num($countsFile)};
	my $maxCount = $countsParse[0];
	my $minCount = $countsParse[1];
	my $maxRange = 200;	# max range allowed in simulations

	my @goodSims; # simulations that passed all 3 tests
	for(my $n = 1; $n <= $simNum*2; $n++){
		my ($rangeTest, $ploidyTest, $overMaxTest) = (0,0,0);
		my $simDir = "$allSimDir/dir_sim_$n/0/";
		my $simMSA = "$simDir/simCounts.txt";
		my $simEvents = "$simDir/simEvents.txt";
		
		# check if maximum count range was reached
		my @chromNumbers;
		open (MSA, "< $simMSA") or die "Can't open simulation MSA file $simMSA\n";
		my @lines = <MSA>;
		foreach my $line (@lines){
			if ($line =~ m/^[^>]/){
				chomp $line;
				push (@chromNumbers, $line);
			}
		}
		my $maxChromNum = max(@chromNumbers);
		my $minChromNum = min(@chromNumbers);
		if ($maxChromNum - $minChromNum < $maxRange){
			$rangeTest = 1;		
		}
		else{
			close MSA;
			next;
		} 
		close MSA;
		
		# check if sim reached max chrom allowed
		open (EVE, "< $simEvents") or die "Can't open simulations events file $simEvents";
		undef $/;
		my $eventsText = <EVE>;
		$eventsText =~ m/Total number of transitions to max chromosome: (\d+)/;
		if ($1 == 0){
			$overMaxTest = 1;
		}
		else{
			close EVE;
			$/ = "\n";
			next;
		}
		close EVE;
		
		# check if all leafs are di/poly ploids
		open (EVE, "< $simEvents") or die "Can't open simulations events file $simEvents";
		$/ = "\n";
		my $eventLine = <EVE>;
		while (defined $eventLine){
			if ($eventLine =~ m/LEAF	GAIN	LOSS	DUPLICATION	DEMI-DUPLICATION/){
				my $leafTotalNum = 0; # total leaf number in simulation
				my $polyPloidNum = 0; # number of di or poly ploid leafs 
				$eventLine = <EVE>;
				while (defined $eventLine){
					$leafTotalNum++;
					$eventLine =~ m/[^\t]+\t\d+\t\d+\t(\d+)\t(\d+)(\t(\d+))?/;
					if ($1 != 0 or $2 != 0 or (defined $4 && $4 != 0)){
						$polyPloidNum++;
					}
					$eventLine = <EVE>;
				}
				if ($polyPloidNum < $leafTotalNum && $polyPloidNum > 0){
					$ploidyTest = 1;
				}
			}
			$eventLine = <EVE>;
		}
		close EVE;
		# decide if sim is OK
		if ($ploidyTest == 1 && $overMaxTest == 1 && $rangeTest == 1){
			push (@goodSims, $n);
		}
	}
	
	# Count good sims
	my $requiredSims = $simNum;
	if (scalar(@goodSims) >= $requiredSims){
		while (scalar(@goodSims) > $requiredSims){
			pop(@goodSims);
		}
		return \@goodSims;	
	}
	else{
		return 0;
	}	
}

# sub prepare fake root frequency file
# creates a root frequency file with lowered largest freq, to be used in fake simulations
sub prepare_fake_root_freq{
	my ($inFreqFile, $fakeRootFreq) = @_;
	
	# parse frequency file to find largest frequency
	my %freqHash;
	open RF , "< $inFreqFile" or die "Can't open $inFreqFile for reading";
	my @frequencies = <RF>;
	foreach my $freq (@frequencies){
		chomp($freq);
		$freq =~ m/F\[(\d+)\]=(.+)/;
		$freqHash{$1} = $2;
	}
	close (RF);
	
	# edit fake root frequency file
	my $largestFreqKey = largest_value_hash({%freqHash});
	my $largestFreq = $freqHash{$largestFreqKey};
	use integer;
	my $keyToSet = $largestFreqKey/2;
	no integer;
	#print LOG "Setting best root frequency to $keyToSet.\n";
	if ($keyToSet < 4){
		print LOG "Best root frequency is too low! Can't reach the required number of required simulations. Process stopped!\n";
		die "Can't reach the required number of good simulations !";
	}
	open RF, "< $inFreqFile" or die "Can't open $inFreqFile for reading";
	my @freqLines = <RF>;
	foreach my $freqLine (@freqLines){
		chomp $freqLine;
		$freqLine =~ m/F\[(\d+)\]/;
		if ($1 == $keyToSet){
			$freqLine =~ s/=.+/=1/;
		}
		else{
			$freqLine =~ s/=.+/=0/;
		}
	}
	close RF;
	open FRF, "> $fakeRootFreq" or die "Can't open $fakeRootFreq for writing";
	print FRF (join("\n", @freqLines));
	close (FRF);
	close LOG;
	return $keyToSet;

}
	
# sub largest value in hash
# Receives a hash (ref) with numeric values, and returns the key with the highest value
sub largest_value_hash {
    my %hash   = %{$_[0]};
    my @keys = keys   %hash;
	my $key = 0;
    my @vals = values %hash;
	my $big = 0;
    for (0 .. $#keys) {
        if ($vals[$_] > $big) {
            $big = $vals[$_];
            $key = $keys[$_];
        }
    }
    $key
}

# sub largest value in hash with median
# Receives a hash (ref) with numeric values and keys, and returns the median
# of the keys with the highest value
sub largest_value_hash_median {
    my %hash   = %{$_[0]};
    my @keys = keys   %hash;
	my $key = 0;
    my @vals = values %hash;
    my @best;
	my $big = 0;
    for (0 .. $#keys) {
        if ($vals[$_] > $big) {
        	$big = $vals[$_];
        	@best = ($keys[$_]); 
        }
        elsif ($vals[$_] == $big){
        	push(@best,$keys[$_]);	
        }       
    }
    my @sortedBest = sort(@best);
    my $median = conf_interval(\@sortedBest,0.5);
    return $median;
}


# sub parse simulation results
# Reads simulated chromosome counts files and prepares them for chromevol inference
sub parse_sim_results{
	### FIELD
	my $facFile = $_[0];
	my $inDir = $_[1];
	my $simListRef = $_[2];
	my $outDir = $_[3];
	
	my @simNumbers = @{$simListRef};
	
	### Run on chosen simulations
	my %refData = %{ &parseFacFile($facFile) };
	foreach my $n (@simNumbers){
		my $inFile = $inDir.'/dir_sim_'.$n.'/0/simCounts.txt';
		### add 'x'
		my %simData = %{ &parseFacFile($inFile) };
		foreach(keys %simData){
			if($refData{$_} =~ m/x/i){
				$simData{$_} = 'x';	# replace simulated count with 'x' if present as 'x' in original chromosome number data
			}
		}
		### print to file
		my $outFile = $outDir.'/sim_'.$n.'.fac';
		open OUT, ">$outFile" or die "ERROR: failed to create $outFile...";
		foreach(sort keys %simData){
			print OUT ">$_\n$simData{$_}\n";
		}
		close OUT;
	}
	### parse FAC file
	sub parseFacFile{
		my $file = $_[0];	# fac file containing chromosome number per species
		my %chrNbr = ();	# store chromosome number by species name
		open FILE, $file or die "ERROR: failed to open $file ...";
		my $name = '';
		my $data = '';
		while(<FILE>){
			chomp;
			my @tmp = split m/\t/, $_;
			if( m/^>/ ){
				$name = substr($_, 1);
			}else{
				$data = $_;
				$chrNbr{$name} = $data;
				$name = '';
				$data = '';
			}
		}
		close FILE;
		return \%chrNbr;
	}
	return 1;
}

# sub prepare inference from simulations
# prepares required files for inferring ploidy from simulated data
sub prepare_infer_from_sim{
	### FIELD
	my $inDirTree = $_[0];
	my $inDirSim = $_[1];
	my $summaryFile = $_[2];
	my $paramFile = $_[3];
	my $goodSimListRef = $_[4];
	my $baseNum = $_[5];
	my $countsFile = $_[6];
	#my $cmdFile = $_[5];
	my $outDir = $_[7];
	my $chromEvolExe = $_[8];
	my @trees = ();  # keep list of trees in Newick format
	
	
	### RUN
	### parse ChromEvol results summary
	# 1) model
	# 2) lnLik
	# 3) AIC
	# 4) dAIC
	open RES, $summaryFile or die "ERROR: failed to open $summaryFile...";
	my $initialParamValues = '';
	while(<RES>){
	        chomp;
	        my @tmp = split m/\t/, $_;
		### get best model
	        if($tmp[3] == 0){
			if($tmp[0] eq 'CONST_RATE'){
				$initialParamValues .= "_gainConstR $tmp[5]\n";
				$initialParamValues .= "_lossConstR $tmp[7]\n";
				$initialParamValues .= "_duplConstR $tmp[9]\n";
			}elsif($tmp[0] eq 'CONST_RATE_DEMI'){
				$initialParamValues .= "_gainConstR $tmp[5]\n";
				$initialParamValues .= "_lossConstR $tmp[7]\n";
				$initialParamValues .= "_duplConstR $tmp[9]\n";
				$initialParamValues .= " _demiPloidyR -2\n";
			}elsif($tmp[0] eq 'CONST_RATE_DEMI_EST'){
				$initialParamValues .= "_gainConstR $tmp[5]\n";
				$initialParamValues .= "_lossConstR $tmp[7]\n";
				$initialParamValues .= "_duplConstR $tmp[9]\n";
				$initialParamValues .= " _demiPloidyR $tmp[10]\n";
			}elsif($tmp[0] eq 'CONST_RATE_NO_DUPL'){
				$initialParamValues .= "_gainConstR $tmp[5]\n";
				$initialParamValues .= "_lossConstR $tmp[7]\n";
			}elsif($tmp[0] eq 'LINEAR_RATE'){
				$initialParamValues .= "_gainConstR $tmp[5]\n";
				$initialParamValues .= "_gainLinearR $tmp[6]\n";
				$initialParamValues .= "_lossConstR $tmp[7]\n";
				$initialParamValues .= "_lossLinearR $tmp[8]\n";
				$initialParamValues .= "_duplConstR $tmp[9]\n";
			}elsif($tmp[0] eq 'LINEAR_RATE_DEMI'){
				$initialParamValues .= "_gainConstR $tmp[5]\n";
				$initialParamValues .= "_gainLinearR $tmp[6]\n";
				$initialParamValues .= "_lossConstR $tmp[7]\n";
				$initialParamValues .= "_lossLinearR $tmp[8]\n";
				$initialParamValues .= "_duplConstR $tmp[9]\n";
				$initialParamValues .= "_demiPloidyR -2\n";
			}elsif($tmp[0] eq 'LINEAR_RATE_DEMI_EST'){
				$initialParamValues .= "_gainConstR $tmp[5]\n";
				$initialParamValues .= "_gainLinearR $tmp[6]\n";
				$initialParamValues .= "_lossConstR $tmp[7]\n";
				$initialParamValues .= "_lossLinearR $tmp[8]\n";
				$initialParamValues .= "_duplConstR $tmp[9]\n";
				$initialParamValues .= "_demiPloidyR $tmp[10]\n";
			}elsif($tmp[0] eq 'LINEAR_RATE_NO_DUPL'){
				$initialParamValues .= "_gainConstR $tmp[5]\n";
				$initialParamValues .= "_gainLinearR $tmp[6]\n";
				$initialParamValues .= "_lossConstR $tmp[7]\n";
				$initialParamValues .= "_lossLinearR $tmp[8]\n";
			}
			elsif($tmp[0] eq 'BASE_NUM'){
				$initialParamValues .= "_gainConstR $tmp[5]\n";
				$initialParamValues .= "_lossConstR $tmp[7]\n";
				$initialParamValues .= "_baseNumberR $tmp[12]\n";
				if ($baseNum != 0){
					$initialParamValues .= "_baseNumber $baseNum\n";
					$initialParamValues .= "_bOptBaseNumber 0\n";
				}
				else{
					# find min count in data - will be used as start point for base number optimization
					my %origCountsHash = %{parse_counts_file($countsFile)};
					my @counts = values(%origCountsHash);
					@counts = grep($_ =~ m/^\d+$/,@counts);	# remove "x" values from sorting
					my @sortCounts = sort {$b <=> $a} @counts;
					my $minInData = $sortCounts[-1];
					my $baseNum = $minInData;
					$baseNum = 4 if $baseNum < 4;
					$initialParamValues .= "_baseNumber $baseNum\n";
					$initialParamValues .= "_bOptBaseNumber 1\n";
					
				}
			}
			elsif($tmp[0] eq 'BASE_NUM_DUPL'){
				$initialParamValues .= "_gainConstR $tmp[5]\n";
				$initialParamValues .= "_lossConstR $tmp[7]\n";
				$initialParamValues .= "_duplConstR $tmp[9]\n";
				$initialParamValues .= "_baseNumberR $tmp[12]\n";
				if ($baseNum != 0){
					$initialParamValues .= "_baseNumber $baseNum\n";
					$initialParamValues .= "_bOptBaseNumber 0\n";
				}
				else{
					# find min count in data - will be used as start point for base number optimization
					my %origCountsHash = %{parse_counts_file($countsFile)};
					my @counts = values(%origCountsHash);
					@counts = grep($_ =~ m/^\d+$/,@counts);	# remove "x" values from sorting
					my @sortCounts = sort {$b <=> $a} @counts;
					my $minInData = $sortCounts[-1];
					my $baseNum = $minInData;
					$baseNum = 4 if $baseNum < 4;
					$initialParamValues .= "_baseNumber $baseNum\n";
					$initialParamValues .= "_bOptBaseNumber 1\n";
					
				}
			}			
	        }
	}
	close RES;
	
	### parse template parameter file for ChromEvol
	open PARAM, $paramFile or die "ERROR: failed to open $paramFile ...";
	undef $/;
	my $paramText = <PARAM>;
	$paramText =~ s/\<PARAMETERS\>/$initialParamValues/;
	$/ = "\n";
	close PARAM;
	
	### print tree and param files
	
	my @goodSims = @{$goodSimListRef};
	
	my @runFilesList;
	foreach my $n (@goodSims){
	        ### print tree to file
	        my $inTreeFile = $inDirTree.'/tree_sim_'.$n;
	        my $inFacFile = $inDirSim.'/sim_'.$n.'.fac';
	        my $outInferDir = $outDir.'/dir_infer_'.$n.'/';
	        ### print param to file
	        my $outParamFile = $outDir.'/param_infer_'.$n;
	        my $newParamText = $paramText;
	        open OUT, ">$outParamFile" or die "ERROR: failed to create $outParamFile...";
	        $newParamText =~ s/\<OUT_DIR\>/$outInferDir/;
	        $newParamText =~ s/\<FAC_FILE\>/$inFacFile/;
	        $newParamText =~ s/\<TREE_FILE\>/$inTreeFile/;
	        print OUT "$newParamText\n";
	        close OUT;
			push(@runFilesList, $outParamFile);
	}
	return \@runFilesList;
}

# sub compute reliability from simulations
# creates a reliability summary from simulations data
sub compute_reliability_from_sim{
	### FIELD
	my $bestModel = $_[0];
	my $resDir = $_[1];
	my $goodSimListRef = $_[2];
	my $inDirInf = $_[3];
	my $inDirSim = $_[4];
	my $thresPP = $_[5];
	my $thresDP = $_[6];
	my $outFile = $_[7];
	my $ploidyCallType = $_[8];
	
	my @goodSimsList = @{$goodSimListRef};

	my $inFileExp = "$resDir/$bestModel/expectations.txt";
		
	my $thresRel = 0.95;
	### get original inference
	my %origExp = %{ &parseExpectationsFile_2thres($inFileExp, $thresDP, $thresPP, $ploidyCallType) };
	### compute total reliability profile
	my %TOTAL_REL = ();
	my $total = 0;
	foreach my $n (@goodSimsList){
		my $INF_FILE = $inDirInf.'/dir_infer_'.$n.'/expectations.txt';
		my $SIM_FILE = $inDirSim.'/dir_sim_'.$n.'/0/simEvents.txt';
		next unless -e $INF_FILE;
		$total++;
		my $INF_EXP = &parseExpectationsFile_2thres($INF_FILE, $thresDP, $thresPP, $ploidyCallType);
		my $SIM_EXP = &parseExpectationsFile($SIM_FILE, 1, $ploidyCallType);
		my %REL_PROFILE = %{ &computeFinalReliabilityProfile($INF_EXP, $SIM_EXP) };
		foreach(keys %REL_PROFILE){
			if(exists $TOTAL_REL{$_}){
				$TOTAL_REL{$_} += $REL_PROFILE{$_};
			}else{
				$TOTAL_REL{$_} = $REL_PROFILE{$_};
			}
		}
	}
	### print profile to file
	open (OUT, "> $outFile") or die "Can't open output file $outFile";
	foreach(sort keys %origExp){
		my $relPercent = $TOTAL_REL{$_}/$total;
		if($relPercent >= $thresRel){
			print OUT "$_\t$origExp{$_}\t$relPercent\n";
		}else{
			print OUT "$_\tNA\t$relPercent\n";
		}
	}
	close OUT;
	return 1;
}

# sub compute reliability from real inference
# creates a reliability summary from real inference data
sub compute_reliability_from_real{
	### FIELD
	my ($inferDir, $inferNum, $thresPP, $thresDP, $outFile,$ploidyCallType) = @_;
	
	# Run on infers, parse expectations file and add to main hash
	my %taxaHash;
	for (my $n = 1; $n <= $inferNum; $n++){
		my $expFile = "$inferDir/dir_$n/expectations.txt";
		my %expHash = %{parseExpectationsFile_2thres($expFile,$thresPP, $thresDP,$ploidyCallType)};
		foreach (keys(%expHash)){
			my @expArr;
			if (defined $taxaHash{$_}){
				@expArr = @{$taxaHash{$_}};
				push (@expArr, $expHash{$_});
			}
			else{
				push (@expArr, $expHash{$_});
			}
			$taxaHash{$_} = [@expArr];
			
		}
	}
	
	# Check reliability for each taxon and print to file
	open (OUT, "> $outFile") or die "Can't open output file $outFile";
	my @taxa = keys(%taxaHash);
	foreach my $taxon (@taxa){
		my @majority = @{majority_in_arr($taxaHash{$taxon})};
		if ($majority[1] >= 0.95){
			print OUT "$taxon\t$majority[0]\t$majority[1]\n"; # if reliable, print ploidy
		}
		else{
			print OUT "$taxon\tNA\t$majority[1]\n"; # if not reliable, print NA
		}
	}
	close(OUT);
	return 1;
}

# sub compute reliability from single topology
# creates a reliability summary from real inference data of single tree topology
sub compute_reliability_from_single{
	### FIELD
	my ($expFile, $thresPP, $thresDP, $outFile,$ploidyCallType) = @_;
	
	# Parse expectations file
	my %expHash = %{parseExpectationsFile_2thres($expFile,$thresPP, $thresDP,$ploidyCallType)};
	
	# Print out
	open (OUT, "> $outFile") or die "Can't open output file $outFile";
	my @taxa = keys(%expHash);
	foreach my $taxon (@taxa){
		print OUT "$taxon\t$expHash{$taxon}\tNA\n";
	}
	close(OUT);
	return 1;
}


# sub majority in array
# receives an array (ref) and returns the most frequent element and it's fraction of the array
sub majority_in_arr{
	my @arr = @{$_[0]};
	my %majHash;
	foreach my $n (@arr){
		if (defined $majHash{$n}){
			$majHash{$n}++;
		}
		else{
			$majHash{$n} = 1;
		}
	}
	
	my $mostFreqKey = largest_value_hash({%majHash});
	my $freq = $majHash{$mostFreqKey};
	my $majority = $freq/scalar(@arr);
	
	my @result = ($mostFreqKey,$majority); # (majority element, fraction)
	
	return [@result];
}

# sub summarize reliability
# creates the final output, based on reliability files previously created
sub summarize_reliability{
	my $relSim = $_[0];
	my $relInf = $_[1];
	my $bestModel = $_[2];
	my $countsFile = $_[3];
	my $outFile = $_[4];
	
	# parse files
	my %relSimHash = %{parse_reliability_file($relSim)};
	my %relInfHash = %{parse_reliability_file($relInf)};
	
	# compare and print to output
	open (OUT, "> $outFile") or die "Can't open output file $outFile";
	print OUT "# Ploidy inference for model $bestModel\n";
	print OUT "# Taxon name\tChromosome count\tPhylogeny robustness\tSimulation reliability\tPloidy inference\n";
	my %countsHash = %{parse_counts_file($countsFile)};
	my @taxa = keys(%relSimHash);
	foreach my $taxon (@taxa){
		my $count = $countsHash{$taxon};
		if ($relInfHash{$taxon}[0] ne "NA" and $relSimHash{$taxon}[0] ne "NA"){
			if ($relInfHash{$taxon}[0] eq "0"){
				print OUT "$taxon\t$count\t$relInfHash{$taxon}[1]\t$relSimHash{$taxon}[1]\t0\n";
			}
			elsif ($relInfHash{$taxon}[0] eq "1"){
				print OUT "$taxon\t$count\t$relInfHash{$taxon}[1]\t$relSimHash{$taxon}[1]\t1\n";
			}
		}
		else{
			print OUT "$taxon\t$count\t$relInfHash{$taxon}[1]\t$relSimHash{$taxon}[1]\tNA\n";
		}
	}
	close OUT;
	return 1;
}

# sub parse reliability file
sub parse_reliability_file{
	my $relFile = $_[0];
	my %resHash;
	open (REL, "< $relFile") or die "Can't open reliability file $relFile";
	my @relLines = <REL>;
	chomp @relLines;
	foreach my $line (@relLines){
		my @splitLine = split(/\t/,$line);
		my @info = ($splitLine[1],$splitLine[2]);
		$resHash{$splitLine[0]} = \@info;
	}
	close REL;
	return {%resHash};
}

# sub get best model
# parses the result summary file to find the best model
sub get_best_model{
	my ($resSum) = @_;
	my $bestModel;
	open RES, $resSum or die "ERROR: failed to open $resSum";
	while(<RES>){
		chomp;
		my @tmp = split m/\t/, $_;
		### get best model
		if($tmp[3] == 0){
			$bestModel = $tmp[0];
		}
	}
	close RES;
	return $bestModel;
}

# sub read input
# reads the control file provided by the user and assigns values to initial parameters
sub read_input{
	my $inputFile = $_[0];
	open (my $in, "< $inputFile") or die "ERROR: Can't read input parameters file";
	my $inGuideTree; my $countsFile; my $inTreesFile; my $infTreeNum = 100;
	my $simNum = 100; my $paramTemplateDir; my $baseNum = 0;
	my $ploidyCallType = "DUPL DEMI BASE"; my $workdir; my $chromEvolExe;
	my $cpusNum = 0; my $duplModelsOnly = 1; my $logFile = "log.txt";
	my $runModels = "ALL"; my $excludeModels = 0; my $plot = 0;
	$baseNum = 0;
	$inTreesFile = "no";
	my @inLines = <$in>;
	chomp @inLines;
	foreach my $line (@inLines){
		if ($line =~ m/^_guideTreeFile (.+)/){
			$inGuideTree = $1;
		}
		elsif ($line =~ m/^_dataFile (.+)/){
			$countsFile = $1;
		}
		elsif ($line =~ m/^_treesFile (.+)/){
			$inTreesFile = $1;
		}
		elsif ($line =~ m/^_topologyNum (.+)/){
			$infTreeNum = $1;
		}
		elsif ($line =~ m/^_simNum (.+)/){
			$simNum = $1;
		}
		elsif ($line =~ m/^_paramTemplates (.+)/){
			$paramTemplateDir = $1;
		}
		elsif ($line =~ m/^_baseNum (.+)/){
			$baseNum = $1;
		}
		elsif ($line =~ m/^_outDir (.+)/){
			$workdir = $1;
			$logFile = "$workdir/log.txt";
		}
		elsif ($line =~ m/^_chromevolExe (.+)/){
			$chromEvolExe = $1;
		}
		elsif ($line =~ m/^_ploidyCallType (.+)/){
			$ploidyCallType = $1;
		}
		elsif ($line =~ m/^_cpusNum (.+)/){
			$cpusNum = $1;
		}
		elsif ($line =~ m/^_ploidyModelsOnly (.+)/){
			$duplModelsOnly = $1;
		}
		elsif ($line =~ m/^_logFile (.+)/){
			$logFile = $1;
		}
		elsif ($line =~ m/^_runModels (.+)/){
			$runModels = $1;
		}
		elsif ($line =~ m/^_excludeModels (.+)/){
			$excludeModels = $1;
		}
		elsif ($line =~ m/^_plot (.+)/){
			$plot = 1;
		}
	}
	my @params = ($inGuideTree,$countsFile,$inTreesFile,$infTreeNum,$simNum,$paramTemplateDir,
		$runModels, $excludeModels, $baseNum,$ploidyCallType,$duplModelsOnly,$workdir,$chromEvolExe,
		$plot,$cpusNum,$logFile);
	foreach (@params) {
		die "ERROR: missing or misspelled parameter $_ in control file"
			unless defined $_;

	}
	return \@params;
}

# sub print params to log
# print the input parameters to the pipeline log file
sub print_params_to_log{
	my @params = @{$_[0]};
	my $logFile = $_[1];
	my %paramsHash = ("_guideTreeFile" => $params[0],"_dataFile" => $params[1],"_treesFile" => $params[2],
	"_topologyNum" => $params[3],"_simNum" => $params[4], "_paramTemplates" => $params[5],
	"_runModels" => $params[6], "_excludeModels" => $params[7], "_baseNum" => $params[8],
	"_ploidyCallType" => $params[9], "_duplModelsOnly" => $params[10],
	"_outDir" => $params[11], "_chromevolExe" => $params[12], "_cpusNum" => $params[13], "_logFile" => $params[14]);
	open (LOG, "> $logFile") or die "ERROR: Can't print parameters to log file '$logFile'\n";
	print LOG "# List of run parameters:\n";	
	foreach (keys(%paramsHash)) {
		print LOG "$_\t$paramsHash{$_}\n";	
	}
	print LOG "\n";
	close LOG;		
}

# sub get ancestral trees
# copies the mlAncestors.tree and posteriorAncestors.tree of best model to the output dir
sub get_ancestral_trees{
	my $workdir = $_[0];
	my $resSumFile = "$workdir/result_sum.txt";
	my $bestModel = get_best_model($resSumFile);
	my $mlAncest = "$workdir/$bestModel/mlAncestors.tree";
	my $posteriorAncest = "$workdir/$bestModel/posteriorAncestors.tree";
	system("cp $mlAncest $workdir/$bestModel\_mlAncestors.tree") == 0
		or die "ERROR: Could not get mlAncestors trees for model $bestModel";
	system("cp $posteriorAncest $workdir/$bestModel\_posteriorAncestors.tree") == 0
		or die "ERROR: Could not get posteriorAncestors trees for model $bestModel";
	return 1;
}

=begin
== This is an old version of the sub, kept in case the new one turns out to make problems ==
# sub create bntpv
# creates the base number transitions probability vector for simulations
sub create_bntpv{
	my ($expFile, $mlAncTree, $baseNum) = @_;	# expectations and  ML Ancestors tree
	
	### parse exp file for nodes with base transition
	my %btnodes;	# all nodes names with base transitions
	open (my $efh, "< $expFile") or die "ERROR: Can't open expectations file $expFile\n";
	my $line = <$efh>;
	while (defined $line){
		chomp $line;
		if ($line =~ m/^#Nodes with transitions in base number/){
			$line = <$efh>;
			while ($line !~ m/Total number/){
				chomp $line;
				$line =~ m/([^:]+): /;
				$btnodes{$1} = 0;
				$line = <$efh>;
			}
			last;
		}
		$line = <$efh>;
	}
	close $efh;
	if (!%btnodes){	# if for some reason there are no base transition events
		return 0;
	}
	### for each node, find the base transition
	my %transitions;	# will contain all transitions found in data (value is # of times transition was found)
	my $treeObj = Bio::TreeIO->new(-file => $mlAncTree, -format => "newick");
	my $tree = $treeObj->next_tree();
	# iterate on all tree nodes
	my $rootnode = $tree->get_root_node;
	foreach my $node ($rootnode->get_all_Descendents()){
		my $nodeid = $node->id();
		$nodeid = $node->bootstrap() if not $nodeid;
		$nodeid =~ m/\[?([^\-]+)-([^\]]+)\]?/;
		my $nodeName = $1; my $nodeCount = $2;
		if (exists $btnodes{$nodeName}){	# if node has base num transition
			my $ancNode = $node->ancestor;	# get parent node's count
			my $ancId = $ancNode->id();
			$ancId = $ancNode->bootstrap() if not $ancId;
			$ancId =~ m/\[?([^\-]+)-(\d+|x)\]?/;
			my $ancCount = $2;
			if ($nodeCount ne "x" && $ancCount ne "x"){	# if both node and parent have counts
				if ($nodeCount =~ m/^\d+$/){	# if normal count
					my $transition = $nodeCount - $ancCount;	# may include gains\losses
					my $realTrans = round($transition/$baseNum)*$baseNum;
					# add transition to hash
					if (exists $transitions{$realTrans}){
						$transitions{$realTrans} += 1;
					}
					else{
						$transitions{$realTrans} = 1;
					}
				}
				else{	# if count is probabilities vector					
					my @parts = split(/_/,$nodeCount);
					foreach my $p (@parts){
						my ($count,$prob) = split(/=/,$p);
						my $transition = $count - $ancCount;	# may include gains\losses
						my $realTrans = round($transition/$baseNum)*$baseNum;
						unless ($realTrans == 0){
							# add transition to hash, weighted by its probability						
							if (exists $transitions{$realTrans}){
								$transitions{$realTrans} += $prob;
							}
							else{
								$transitions{$realTrans} = $prob;
							}
						}
					}
				}
			}
			else{
				$btnodes{$nodeName} = "NA";
			}
		}

	}
	### calculate probabilities
	my @allTransitions = values(%transitions);
	my $totalTrans = eval(join('+',@allTransitions));
	# remove rare transitions (less then 0.05)
	foreach my $t (keys(%transitions)){
		if ($transitions{$t}/$totalTrans < 0.05){
			delete $transitions{$t};
		}
	}
	# create final probs hash
	my %probsHash;
	@allTransitions = values(%transitions);
	$totalTrans = eval(join('+',@allTransitions));
	foreach my $t (keys(%transitions)){
		my $finalProb = sprintf("%.2f",$transitions{$t}/$totalTrans);	# limited to 2 decimal digits
		$probsHash{$t} = $finalProb;
	}
	# make sure probs sum to 1
	my @allProbs = values(%probsHash);
	my $probSum = eval(join('+',@allProbs));
	my $diff = 1 - $probSum;
	if ($diff != 0){
		# complete or substract to 1 on a random prob
		my @transitions = keys(%probsHash);
		my $randTrans = $transitions[0];
		$probsHash{$randTrans} += $diff;
	}

	# create probs vector
	my $vector = "";
	foreach my $p (keys(%probsHash)){
		$vector .= $p."=".$probsHash{$p}."_"
	}
	$vector = substr($vector,0,-1);
	return $vector;	
}
=cut

# sub create bntpv
# creates the base number transitions probability vector for simulations
sub create_bntpv{
	my ($expFile, $mlAncTree, $baseNum) = @_;	# expectations and  ML Ancestors tree
	
	### parse exp file for nodes with base transition > 0.1
	my %btnodes;	# all nodes names with base transitions
	open (my $efh, "< $expFile") or die "ERROR: Can't open expectations file $expFile\n";
	my $line = <$efh>;
	while (defined $line){
		chomp $line;
		if ($line =~ m/^#ALL EVENTS EXPECTATIONS PER NODE/){
			$line = <$efh>;
			$line = <$efh>;
			while ($line !~ m/#\++/){
				chomp $line;
				my @transitions = split(/\t/,$line);
				my $node = $transitions[0];
				my $weight = $transitions[5];
				if ($weight > 0.1){
					$weight = 1 if $weight > 1;
					$btnodes{$node} = $weight;	
				}				
				$line = <$efh>;
			}
			last;
		}
		$line = <$efh>;
	}
	close $efh;
	if (!%btnodes){	# if for some reason there are no base transition events
		return "$baseNum\=1.00";
	}
	
	### for each node, find the base transition
	my %transitions;	# will contain all transitions found in data (value is # of times transition was found)
	my $treeObj = Bio::TreeIO->new(-file => $mlAncTree, -format => "newick");
	my $tree = $treeObj->next_tree();
	# iterate on all tree nodes
	my $rootnode = $tree->get_root_node;
	foreach my $node ($rootnode->get_all_Descendents()){
		my $nodeid = $node->id();
		$nodeid = $node->bootstrap() if not $nodeid;
		$nodeid =~ m/\[?([^\-]+)-([^\]]+)\]?/;
		my $nodeName = $1; my $nodeCount = $2;
		if (exists $btnodes{$nodeName}){	# if node has base num transition
			my $ancNode = $node->ancestor;	# get parent node's count
			my $ancId = $ancNode->id();
			$ancId = $ancNode->bootstrap() if not $ancId;
			$ancId =~ m/\[?([^\-]+)-(\d+|x)\]?/;
			my $ancCount = $2;
			if (($nodeCount ne "x" and $nodeCount ne "X") && ($ancCount ne "x" and $ancCount ne "X")){	# if both node and parent have counts
				if ($nodeCount =~ m/^\d+$/){	# if normal count
					my $transition = $nodeCount - $ancCount;	# may include gains\losses
					my $realTrans = round($transition/$baseNum)*$baseNum;
					# add transition to hash
					if ($realTrans > 0){
						if (exists $transitions{$realTrans}){
							$transitions{$realTrans} += $btnodes{$nodeName};
						}
						else{
							$transitions{$realTrans} = $btnodes{$nodeName};
						}
					}
				}
				elsif ($nodeCount eq 'x' or $nodeCount eq 'X'){
					next;
				}
				else{	# if count is probabilities vector					
					my @parts = split(/_/,$nodeCount);
					foreach my $p (@parts){
						my ($count,$prob) = split(/=/,$p);
						my $transition = $count - $ancCount;	# may include gains\losses
						my $realTrans = round($transition/$baseNum)*$baseNum;
						if ($realTrans > 0){					
							if (exists $transitions{$realTrans}){
								$transitions{$realTrans} += $prob*$btnodes{$nodeName};
							}
							else{
								$transitions{$realTrans} = $prob*$btnodes{$nodeName};
							}
						}
					}
				}
			}
			else{
				$btnodes{$nodeName} = "NA";
			}
		}

	}
	# in case no transitions were found, give transition by base number prob of 1
	if (not %transitions){
		$transitions{$baseNum} = 1.00; 
	}
	
	### calculate probabilities
	my @allTransitions = values(%transitions);
	my $totalTrans = eval(join('+',@allTransitions));
	# remove rare transitions (less then 0.05)
	foreach my $t (keys(%transitions)){
		if ($transitions{$t}/$totalTrans < 0.05){
			delete $transitions{$t};
		}
	}
	# create final probs hash
	my %probsHash;
	@allTransitions = values(%transitions);
	$totalTrans = eval(join('+',@allTransitions));
	foreach my $t (keys(%transitions)){
		my $finalProb = sprintf("%.2f",$transitions{$t}/$totalTrans);	# limited to 2 decimal digits
		$probsHash{$t} = $finalProb;
	}
	# make sure probs sum to 1
	my @allProbs = values(%probsHash);
	my $probSum = eval(join('+',@allProbs));
	my $diff = 1 - $probSum;
	if ($diff != 0){
		# complete or substract to 1 on a random prob
		my @transitions = keys(%probsHash);
		my $randTrans = $transitions[0];
		$probsHash{$randTrans} += $diff;
	}

	# create probs vector
	my $vector = "";
	foreach my $p (keys(%probsHash)){
		$vector .= $p."=".$probsHash{$p}."_"
	}
	$vector = substr($vector,0,-1);
	return $vector;	
}

# sub round
# receives a float number and returns rounded number
sub round{
	my $n = $_[0];
	my $roundup = int($n+0.5);
	my $rounddown = int($n-0.5);
	if ($n - $rounddown < $roundup - $n){
		return $rounddown;
	}
	else{
		return $roundup;
	}
}


# sub determine thresholds
# Based on simulated data, this unction calculates the optimal thresholds
# for calling diploids and polyploids (two different thresholds)
sub determine_thresholds{
	my $goodSimListRef = $_[0];
	my $inDirInf = $_[1];
	my $inDirSim = $_[2];
	my $outFilePP = $_[3];
	my $outFileDP = $_[4];
	my $ploidyCallType = $_[5];
	
	## determine PP threshold
	my $thresPP = determine_PP_threshold($goodSimListRef,$inDirInf,$inDirSim,
		$outFilePP, $ploidyCallType);
	## determine DP threshold
	my $thresDP = determine_DP_threshold($goodSimListRef,$inDirInf,$inDirSim,
		$outFileDP, $ploidyCallType, $thresPP);	# PP threshold given as argument
	
	my @thresholds = ($thresPP,$thresDP);
	return \@thresholds;
}

# sub MCC
# calculates Mathew's Correlation Coefficient
sub MCC{
	my ($TP, $TN, $FP, $FN) = @_;
	my $MCC = 0;
	# if MCC denominator is zero, set MCC value to zero
	if (($TP+$FP)*($TP+$FN)*($TN+$FP)*($TN+$FN) != 0){
		$MCC = ($TP*$TN - $FP*$FN)/(($TP+$FP)*($TP+$FN)*($TN+$FP)*($TN+$FN))**0.5;
	}	
	return $MCC;
}

# sub determine Polyploid threshold
# determines the optimal threshold for calling polyploids
sub determine_PP_threshold{
	my $goodSimListRef = $_[0];
	my $inDirInf = $_[1];
	my $inDirSim = $_[2];
	my $outFile = $_[3];
	my $ploidyCallType = $_[4];
	my @goodSimList = @{$goodSimListRef};
	open (my $thresfh, "> $outFile") or die "Can't open output thresholds file $outFile";
	
	my %thresMCC;
	for(my $thresP = 0; $thresP < 1.01; $thresP+=0.01){
		my $thresResultsRef = test_threshold($goodSimListRef,$inDirInf,$inDirSim,
			$outFile,$ploidyCallType, $thresP,1);
		my %thresResults = %{$thresResultsRef};
		my ($thresTP, $thresTN, $thresFP, $thresFN) = 
			($thresResults{"TP"}, $thresResults{"TN"}, $thresResults{"FP"}, $thresResults{"FN"});
		my $mcc = MCC($thresTP, $thresTN, $thresFP, $thresFN);
		$thresMCC{$thresP} = $mcc;
		print $thresfh "$thresP\t$mcc\n";
	}
	my $bestThresP = largest_value_hash_median(\%thresMCC);
	print $thresfh "#Best threshold: $bestThresP\n";
	close $thresfh;
	return $bestThresP;
}

# sub determine Diploid threshold
# determines the optimal threshold for calling diploids
sub determine_DP_threshold{
	my $goodSimListRef = $_[0];
	my $inDirInf = $_[1];
	my $inDirSim = $_[2];
	my $outFile = $_[3];
	my $ploidyCallType = $_[4];
	my $thresP = $_[5];	# PP threshold given as upper limit or DP threshold
	my @goodSimList = @{$goodSimListRef};
	open (my $thresfh, "> $outFile") or die "Can't open output thresholds file $outFile";
	
	my %thresMCC;
	for(my $thresD = 0; $thresD < $thresP; $thresD+=0.01){
		my $thresResultsRef = test_threshold($goodSimListRef,$inDirInf,$inDirSim,
			$outFile,$ploidyCallType, $thresD,0);
		my %thresResults = %{$thresResultsRef};
		my ($thresTP, $thresTN, $thresFP, $thresFN) = 
			($thresResults{"TP"}, $thresResults{"TN"}, $thresResults{"FP"}, $thresResults{"FN"});
		
		my $mcc = MCC($thresTP, $thresTN, $thresFP, $thresFN);
		$thresMCC{$thresD} = $mcc;
		print $thresfh "$thresD\t$mcc\n";
	}
	my $bestThresD = largest_value_hash_median(\%thresMCC);
	print $thresfh "#Best threshold: $bestThresD\n";
	close $thresfh;
	return $bestThresD;
}

# sub test threshold
# given a threshold, the function returns a hash with the total
# numbers of TP,TN,FP,FN in all simulations
sub test_threshold{
	my $goodSimListRef = $_[0];
	my $inDirInf = $_[1];
	my $inDirSim = $_[2];
	my $outFile = $_[3];
	my $ploidyCallType = $_[4];
	my $thres = $_[5];	# threshold value to be tested
	my $type = $_[6];	# 0 ('under') for DP, 1 ('above') for PP
	my @goodSimList = @{$goodSimListRef};
	
	my %testResults = ("TP"=>0, "TN"=>0, "FP"=>0, "FN"=>0);	
	### iterate through all simulations
	foreach my $n (@goodSimList){
		### make file name
		my $INF_FILE = $inDirInf.'/dir_infer_'.$n.'/expectations.txt';
		my $SIM_FILE = $inDirSim.'/dir_sim_'.$n.'/0/simEvents.txt';
		next unless -e $INF_FILE;
		my $INF_EXP = &parseExpectationsFile($INF_FILE, $thres, $ploidyCallType);		
		my $SIM_EXP = &parseExpectationsFile($SIM_FILE, 1, $ploidyCallType);
		my %result = %{ &computeReliabilityProfile($INF_EXP, $SIM_EXP, $type) };
		### get reliability profile
		$testResults{"TP"} += $result{"TP"};
		$testResults{"TN"} += $result{"TN"};
		$testResults{"FP"} += $result{"FP"};
		$testResults{"FN"} += $result{"FN"};
	}
			
	return \%testResults;
	
}

# sub parse expectations file
#returns expected number of events from root to leaf.
# Can either work on an inferred expectations file or simulation events file

#LEAF    GAIN    LOSS    DUPLICATION     DEMI-DUPLICATION	(BASE_NUMBER)
#Onychostoma_gerlachi    0       1.01103 0       0
#Tor_macrolepis  0       0.0509317       1.84544 0.000364298
#Tor_putitora    0       0.0111298       1.99996 6.806e-07
sub parseExpectationsFile{
	my $file = $_[0];
	my $thres = $_[1];	# cutoff threshold
	my $ploidyCallType = $_[2];
	my %EXP = ();		# store expectations indexed by species name
	open FILE, $file or die "ERROR: failed to open $file ...";
	my $parse = 0;	# parse only when >= 2
	while(<FILE>){
		chomp;
		if( m/\#Expected number of events from root to leaf/ || m/LEAF\s+GAIN\s+LOSS\s+DUPLICATION\s+DEMI-DUPLICATION/ ){
			$parse++;
			next;
		}
		if( $parse >= 2 ){
			my @tmp = split m/\t/, $_;
			my $duplEvents = $tmp[3];
			my $demiduplEvents = $tmp[4];
			my $basenumEvents = $tmp[5];
			# calculate number of events
			my $events = 0;
			my @requiredEvents = split(/\s/, $ploidyCallType);
			foreach my $callType (@requiredEvents){
				$events += $duplEvents if (defined $duplEvents && $callType eq "DUPL");
				$events += $demiduplEvents if (defined $demiduplEvents && $callType eq "DEMI");
				$events += $basenumEvents if (defined $basenumEvents && $callType eq "BASE");
			}
			if( $events >= $thres ){
				$EXP{$tmp[0]} = 1;
			}elsif( $events < $thres ){
				$EXP{$tmp[0]} = 0;
			}
		}
	}
	delete $EXP{''};
	close FILE;
	return \%EXP;
}

# sub compute reliability profile
# compares simulated ploidy against ploidy inferred from simulations
sub computeReliabilityProfile{
	my %INF = %{ $_[0] };
	my %SIM = %{ $_[1] };
	my $type = $_[2];
	my %REL = ("TP"=>0, "TN"=>0, "FP"=>0, "FN"=>0);
	if ($type == 1){	# PP threshold
		foreach(keys %INF){
			if ($INF{$_} == $SIM{$_} && $INF{$_} == 1){
				$REL{"TP"} += 1;
			}
			elsif($INF{$_} == $SIM{$_} && $INF{$_} == 0){
				$REL{"TN"} += 1;
			}
			elsif($INF{$_} != $SIM{$_} && $INF{$_} == 1){
				$REL{"FP"} += 1;
			}
			elsif($INF{$_} != $SIM{$_} && $INF{$_} == 0){
				$REL{"FN"} += 1;
			}
		}
	}
	elsif ($type == 0){	# DP threshold
		foreach(keys %INF){
			if ($INF{$_} == $SIM{$_} && $INF{$_} == 0){
				$REL{"TP"} += 1;
			}
			elsif($INF{$_} == $SIM{$_} && $INF{$_} == 1){
				$REL{"TN"} += 1;
			}
			elsif($INF{$_} != $SIM{$_} && $INF{$_} == 0){
				$REL{"FP"} += 1;
			}
			elsif($INF{$_} != $SIM{$_} && $INF{$_} == 1){
				$REL{"FN"} += 1;
			}
		}	
	}
	delete $REL{''};
	return \%REL;
}


# sub get params distribution
# creates an output of rate parameters distribution and confidence interval
sub get_params_distribution{
	my ($inferDir, $topologyNum, $outFile) = @_;
	my (@gainConst, @lossConst, @dupl, @demi, @baseNum, @gainLinear, @lossLinear);
	## collect rates
	for (my $i=1; $i<=$topologyNum; $i++){
		my $resFile = "$inferDir/dir_$i/chromEvol.res";
		open (my $rfh, "< $resFile") or die "ERROR: Can't open $resFile\n";
		my @resLines = <$rfh>;
		chomp @resLines;
		close $rfh;
		foreach my $line (@resLines){
			if ($line =~ m/LOSS_CONST\t(.+)/){
				push(@lossConst,$1);
			}
			elsif ($line =~ m/GAIN_CONST\t(.+)/){
				push(@gainConst,$1);
			}
			elsif ($line =~ m/DUPL\t(.+)/){
				push(@dupl,$1);
			}
			elsif ($line =~ m/BASE_NUMBER_R\t(.+)/){
				push(@baseNum,$1);
			}
			elsif ($line =~ m/HALF_DUPL\t(.+)/){
				push(@demi,$1);
			}
			elsif ($line =~ m/LOSS_LINEAR\t(.+)/){
				push(@lossLinear,$1);
			}
			elsif ($line =~ m/GAIN_LINEAR\t(.+)/){
				push(@gainLinear,$1);
			}
		}
	}
	
	## summarize distribution
	open (my $ofh, "> $outFile") or die "ERROR: CAn't write to parameter distribution file $outFile\n";
	print $ofh "Parameter\tmedian(90% interval)\n";
	print $ofh "GAIN_CONST_RATE\t".summarize_dist(\@gainConst)."\n" if (scalar(@gainConst)== 100);
	print $ofh "LOSS_CONST_RATE\t".summarize_dist(\@lossConst)."\n" if (scalar(@lossConst)== 100);
	print $ofh "DUPL_RATE\t".summarize_dist(\@dupl)."\n" if (scalar(@dupl)== 100);
	print $ofh "DEMI_DUPL_RATE\t".summarize_dist(\@demi)."\n" if (scalar(@demi)== 100);
	print $ofh "BASE_TRANSITION_RATE\t".summarize_dist(\@baseNum)."\n" if (scalar(@baseNum)== 100);
	print $ofh "GAIN_LINEAR_RATE\t".summarize_dist(\@gainLinear)."\n" if (scalar(@gainLinear)== 100);
	print $ofh "LOSS_LINEAR_RATE\t".summarize_dist(\@lossLinear)."\n" if (scalar(@lossLinear)== 100);
	close $ofh;
	return 1;	

	
	sub summarize_dist{
		my @paramArr = @{$_[0]};
		my @sortedParamArr = sort { $a <=> $b } @paramArr;
		my $median = conf_interval(\@sortedParamArr, 0.5);
		my $low = conf_interval(\@sortedParamArr, 0.05);
		my $high = conf_interval(\@sortedParamArr, 0.95);
		return "$median ($low\-$high)";
	}
}

# sub confidence interval
# recieves an array (reference) and an interval (between 0 and 1) and returns the corresponding percentile
sub conf_interval{
	my @sortedArr = @{$_[0]};
	my $interval = $_[1];	# 0 to 1
	my $n = scalar(@sortedArr);
	my $index = round($n*$interval);
	$index-- if $index > scalar(@sortedArr) - 1;	# prevent indexes out of range
	my $bound = $sortedArr[$index];
	return $bound;
}

# sub prepare simulations ### !!!!!!! not working properly. FIX!
# prepares parameter, trees and commands files for simulations
sub prepare_simulations{
	### FILES AND DIR
	my ($infDir, $freqFile, $paramTemplate, $nbrSims, $maxChrNum,
		$chromNumRange, $simDir, $chromEvolExe, $fake, $single) = @_;
	
	system("rm -rf $simDir/*") if $fake;	# if fake simulations, empty simulations dir
	### FIELD
	my @trees = ();	 # keep list of trees in Newick format
	my $minChrNum = 1;	# CHANGE if needed
	my @simParamFilesList;
	
	my $count = 1;
	my $n = 1;
	my $singleRes; my $singleExp; my $singleMlAnc;
	if ($single){	# if only single tree, take parameters from initial run
		my $workdir = dirname($infDir);
		my $resSum = "$workdir/result_sum.txt";
		my $bestModel = get_best_model($resSum);
		$singleRes = "$workdir/$bestModel/chromEvol.res";
		$singleExp = "$workdir/$bestModel/expectations.txt";
		my $singleMlAnc = "$workdir/$bestModel/mlAncestors.tree";
	}

	while ($count <= $nbrSims*2){
		my $simNum = $count;
		my $tree = "$infDir/tree_$n";
		$tree = $single if $single;	# if only single tree available, use it
		my $outTreeFile = "$simDir/tree_sim_$simNum";
		if (-e $tree){
			copy($tree, $outTreeFile);
			# get param values
			my $totalTreeLength;
			my $paramValues = "";
			my $resFile = "$infDir/dir_$n/chromEvol.res";
			$resFile = $singleRes if $single;	# if only single tree, take parameters from initial run
			open (my $rfh, "< $resFile") or die "ERROR: Can't open $resFile\n";
			my @resLines = <$rfh>;
			chomp @resLines;
			close $rfh;
			my $baseNum;
			foreach my $line (@resLines){
				if ($line =~ m/LOSS_CONST\t(.+)/){
					$paramValues .= "_lossConstR $1\n";
				}
				elsif ($line =~ m/GAIN_CONST\t(.+)/){
					$paramValues .= "_gainConstR $1\n";
				}
				elsif ($line =~ m/DUPL\t(.+)/){
					$paramValues .= "_duplConstR $1\n";
				}
				elsif ($line =~ m/BASE_NUMBER_R\t(.+)/){
					$paramValues .= "_baseNumberR $1\n";
				}
				elsif ($line =~ m/BASE_NUMBER\t(.+)/){
					$paramValues .= "_baseNumber $1\n";
					$baseNum = $1;
				}		
				elsif ($line =~ m/HALF_DUPL\t(.+)/){
					$paramValues .= "_demiPloidyR $1\n";
				}
				elsif ($line =~ m/LOSS_LINEAR\t(.+)/){
					$paramValues .= "_lossLinearR $1\n";
				}
				elsif ($line =~ m/GAIN_LINEAR\t(.+)/){
					$paramValues .= "_gainLinearR $1";
				}
				elsif ($line =~ m/#total tree length = (.+)/){
					$totalTreeLength = $1;
				}				
			}
			close $rfh;
			if ($baseNum){	# if base number model, get bntpv
				my $expFile = "$infDir/dir_$n/expectations.txt";
				$expFile = $singleExp if $single;
				my $mlAncTree = "$infDir/dir_$n/mlAncestors.tree";
				$mlAncTree = $singleMlAnc if $single;
				my $bntpv = create_bntpv($expFile, $mlAncTree, $baseNum);
				die "ERROR: Can't acquire base number transitions probailities vector for simulations!\n" unless $bntpv;
				$paramValues .= "_baseTransitionProbs $bntpv";
			}
			
			### parse parameter file for ChromEvol simulation
			open PARAM, $paramTemplate or die "ERROR: failed to open $$paramTemplate...";
			undef $/;
			$maxChrNum = 200 if ($fake && $maxChrNum < 200);	# if creating fake simulations
			my $paramText = <PARAM>;
			$paramText =~ s/\<FREQ_FILE\>/$freqFile/;
			$paramText =~ s/\<MIN_CHR_NUM\>/$minChrNum/;
			$paramText =~ s/\<MAX_CHR_NUM\>/$maxChrNum/;
			$paramText =~ s/\<PARAMETERS\>/$paramValues/;
			$paramText =~ s/\<TREE_LENGTH\>/$totalTreeLength/;
			$paramText =~ s/\<MAX_CHR_NUM_SIM\>/$maxChrNum/;
			
			
			$/ = "\n";
			close PARAM;
			
			### make directory to contain simulations results
			my $outSimDir = $simDir.'/dir_sim_'.$simNum;
			### print param to file
			my $outParamFile = $simDir.'/param_sim_'.$simNum;
			my $newParamText = $paramText;
			open OUT, ">$outParamFile" or die "ERROR: failed to create $outParamFile...";
			$newParamText =~ s/\<OUT_DIR\>/$outSimDir/;
			$newParamText =~ s/\<TREE_FILE\>/$outTreeFile/;
			print OUT $newParamText."\n";
			close OUT;
			### add param file to list
			push(@simParamFilesList, $outParamFile);
			
			$n++;
			$count++;	
		}
		else{
			$n = 1;
		}	
	}
	
	return \@simParamFilesList;
}

# sub parse expectations file with two thresholds
# similar to parse expectations file, but considers a different threshold for DP and PP
sub parseExpectationsFile_2thres{
	my $file = $_[0];
	my $thresPP = $_[1];
	my $thresDP = $_[2];
	my $ploidyCallType = $_[3];
	my %EXP = ();		# store expectations indexed by species name
	open FILE, $file or die "ERROR: failed to open $file ...";
	my $parse = 0;	# parse only when >= 2
	while(<FILE>){
		chomp;
		if( m/\#Expected number of events from root to leaf/ || m/LEAF\s+GAIN\s+LOSS\s+DUPLICATION\s+DEMI-DUPLICATION/ ){
			$parse++;
			next;
		}
		if( $parse >= 2 ){
			my @tmp = split m/\t/, $_;
			my $duplEvents = $tmp[3];
			my $demiduplEvents = $tmp[4];
			my $basenumEvents = $tmp[5];
			# calculate number of events
			my $events = 0;
			my @requiredEvents = split(/\s/, $ploidyCallType);
			foreach my $callType (@requiredEvents){
				$events += $duplEvents if (defined $duplEvents && $callType eq "DUPL");
				$events += $demiduplEvents if (defined $demiduplEvents && $callType eq "DEMI");
				$events += $basenumEvents if (defined $basenumEvents && $callType eq "BASE");
			}
			if( $events >= $thresPP ){
				$EXP{$tmp[0]} = 1;			
			}
			elsif ($events > $thresDP && $events < $thresPP){
				$EXP{$tmp[0]} = "NA";
			}
			elsif( $events <= $thresDP ){
				$EXP{$tmp[0]} = 0;
			}
		}
	}
	delete $EXP{''};
	close FILE;
	return \%EXP;
}

# sub compute final reliability profile
# compares simulated ploidy against ploidy inferred from simulations
sub computeFinalReliabilityProfile{
	my %INF = %{ $_[0] };
	my %SIM = %{ $_[1] };
	my %REL;
	foreach(keys %INF){
		if ($INF{$_} == $SIM{$_}){
			if (exists $REL{$_}){
				$REL{$_} += 1;
			}
			else{
				$REL{$_} = 1;
			}
		}
		else{
			if ( not exists $REL{$_}){
				$REL{$_} += 0;
			}			
		}
	}
	
	delete $REL{''};
	return \%REL;
}

# sub choose trees
# choses a set of trees out of the given trees. Returns a ref to list of trees numbers
sub choose_trees{
	my ($treesFile, $choiceSize) = @_;
	open (my $tfh, "< $treesFile") or die "ERROR: Can't open input trees file\n";
	my @trees = <$tfh>;
	close $tfh;
	my $treesInFile = scalar(@trees);
	if ($treesInFile < $choiceSize){
		warn "ERROR: Number of required trees is larger than number of trees in input file\n";
		return 0;
	}
	my @chosen;	
	my @shuffled_indexes = shuffle(0..$treesInFile - 1);
	# pick trees at random
	my $tr = 0;
	while (scalar(@chosen) < $choiceSize){
		if ($trees[$shuffled_indexes[$tr]]){	# make sure line contains a tree
			push(@chosen, $shuffled_indexes[$tr]);
		}
		$tr++;		
	}		

	return \@chosen;
}


# sub print good simulations list
# prints the chosen simulations numbers to a file
sub print_good_sims{
	my ($goodSimListRef, $workdir) = @_;
	my @goodSimList = @{$goodSimListRef};
	my $outfile = "$workdir/chosen_sims";
	open (my $fh, "> $outfile") or die "ERROR: Can't write chosen simulations to file $outfile\n";
	print $fh join("\n",@goodSimList);
	close $fh;
}


# sub check runs
# checks if all intended initial runs were completed successfully
sub check_runs{
	my ($runFilesRef, $workdir) = @_;
	my @runFiles = @{$runFilesRef};
	my $intendedRuns = scalar(@runFiles);
	my $completedCount = 0;
	foreach my $file (@runFiles){
		my ($modelDir) = ($file =~ m/param_(.+)\.txt/);
		if (-e "$workdir/$modelDir/chromEvol.res"){
			$completedCount++;
		}
	}
	if ($completedCount == $intendedRuns){
		return 1;
	}
	else{
		return 0;
	}
}

1;