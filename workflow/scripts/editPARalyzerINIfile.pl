( $iniFile, $outputDirectory, $inputFile, $tbitFile ) = @ARGV;

open(USRFL, "<", $iniFile) or die "Can't open $iniFile\n";
while(<USRFL>) {
	if( not /^[BOWTIE_FILE|OUTPUT]/ && length($_) > 0 ) {
		print $_;
	}
}
close USRFL;

print "\n";
print "GENOME_2BIT_FILE=$tbitFile\n";
print "SAM_FILE=$inputFile\n\n";
print "OUTPUT_DISTRIBUTIONS_FILE=$outputDirectory\.distribution\n";
print "OUTPUT_GROUPS_FILE=$outputDirectory\.groups\n";
print "OUTPUT_CLUSTERS_FILE=$outputDirectory\.clusters\n";
print "OUTPUT_READS_FILE=$outputDirectory\_PARalyzer_Utilized.sam\n\n";