 
// Sliding dynamic window binning.
// Remove duplicates
// 1% category

// import modules
import std.getopt, std.math, std.stdio, std.regex, std.algorithm, std.range, std.random, std.file, std.string, std.conv;
import bio.bam.reader, bio.bam.pileup;


// initialize global variables
bool help;
string patientFolder;
string fatherFolder;
string motherFolder;
int MQ = 30;
ulong reads_per_bin = 15;
double distance_factor = 3;
double patient_percentage = 65;
double patient_count_perc = 20;

int[string][string] qc_fail;
double patient_count;
double father_count;
double mother_count;
string patientName;
string fatherName;
string motherName;

// function to show usage/help
void printUsage() {
    stderr.writeln("Usage: switch_point_caller_trio [options] <input.bam> ");
    stderr.writeln();
    stderr.writeln("Options: -h, --help");
    stderr.writeln("                    show usage");        
    stderr.writeln("         -p, --patient (required)");
    stderr.writeln("                    folder with the patient bam files");    
    stderr.writeln("         -f, --father (required)");
    stderr.writeln("                    folder with the father bam files");    
    stderr.writeln("         -m, --mother (required)");
    stderr.writeln("                    folder with the mother bam files");        
    stderr.writeln("         -r, --reads");
    stderr.writeln("                    number of reads per bin (DEFAULT: 15)");
    stderr.writeln("         -d, --distance");
    stderr.writeln("                    distance factor (DEFAULT: 3)");
    stderr.writeln("         -x, --percentage");
    stderr.writeln("                    minimal percentage of patient cells of all cells supporting the switchpoint(DEFAULT: 65%)");
    stderr.writeln("         -c, --count");
    stderr.writeln("                    minimal percentage of 'passed' cells in the patient supporting the switchpoint (DEFAULT: 20%)");
    stderr.writeln("         -q, --MQ");
    stderr.writeln("                    use only reads with high mapping quality (DEFAULT: 30)");
}

int main(string[] args) {
    getopt(args, std.getopt.config.caseSensitive, std.getopt.config.passThrough,
      "help|h", &help,
      "reads|r", &reads_per_bin,    
      "distance|d", &distance_factor,
      "percentage|x", &patient_percentage,
      "count|c", &patient_count_perc,
      "MQ|q", &MQ,
      "patient|p", &patientFolder,
      "father|f", &fatherFolder,
      "mother|m", &motherFolder,    
    );
  
  // die and show usage if help or patient is empty or father is empty or mother is empty
  if (help || patientFolder == null || fatherFolder == null || motherFolder == null) {
    printUsage();
    return 0;
  }
  
  // get name of the sample
  patientName = get_output_name(patientFolder);
  fatherName = get_output_name(fatherFolder);
  motherName = get_output_name(motherFolder);
  
  // get median distributions
  double[] bin_size_median_distribution_patient = bam_looping(patientFolder);
  double[] bin_size_median_distribution_father = bam_looping(fatherFolder);
  double[] bin_size_median_distribution_mother = bam_looping(motherFolder);
  
  double[] bin_size_median_distribution_trio = bin_size_median_distribution_patient ~ bin_size_median_distribution_father ~ bin_size_median_distribution_mother; // merge median distributions
  double bin_size_median_trio_mean = mean(bin_size_median_distribution_trio); // calculate the mean of the median distribution
  double merge_distance = bin_size_median_trio_mean/distance_factor; // calculate the maximal merge distance between two switchpoints
  merge_switch_points(patientFolder, fatherFolder, motherFolder, merge_distance); // call denovo switch points clusters
  
  return 1;
}

// function to loop through the folder and search for *.bam files.
double[] bam_looping( string folder ) {
  // initialize variables
  double[] total_reads_distribution;
  double[] FR_bins_ratio_distribution;
  double[] bin_size_median_distribution;
  double[] number_of_bins_distribution;
  
  auto bamFiles = array(filter!`endsWith(a.name,".bam")`(dirEntries(folder,SpanMode.depth))); // find all *.bam files in the folder
    
  foreach( bamFile ; bamFiles ) { // loop through the bam files
    auto output = process_bam(bamFile); // process a single bam file
      
    // store output in different variable
    total_reads_distribution ~= output[0];
    FR_bins_ratio_distribution ~= output[1];
    bin_size_median_distribution ~= output[2];
    number_of_bins_distribution ~= output[3];
  }  
  bin_size_median_distribution = qc_analysis(bamFiles, total_reads_distribution, FR_bins_ratio_distribution, bin_size_median_distribution, number_of_bins_distribution, folder); // preform QC analysis
  return(bin_size_median_distribution);  // return median distribution
}

// function to process a single bam file
Tuple!(double, double, double, double) process_bam(string file) {
  // initialize variables
  double total_reads = 0;
  double number_of_bins = 0;
  double[] bin_size_distribution;
  ulong FR_bins;
  ulong total_bins;    
  
  stdout.writeln("Busy with "~file~"..."); // show progress  
  BamReader bam = new BamReader(file); // open bam file with the bamreads module
    
  // open output files
  File binsFile = File(file~".bins","w");
  File regionsFile = File(file~".regions","w");
  File pointsFile = File(file~".points","w");
  
  // write headers to the output files
  binsFile.writeln("#CHR\tSTART\tEND\tSIZE\tF\tR\tCAT");
  regionsFile.writeln("#CHR\tSTART\tEND\tSIZE\tCAT");
  pointsFile.writeln("#CHR\tPOS");
  
  foreach( chromosome ;  bam.reference_sequences() ) { // loop through each chromosome in a single bam file
    if (startsWith(chromosome.name, "GL")) { continue; } // skip reads on chromosome contigs
    auto output_determine_bins = determine_bins(bam, chromosome); // determine bins based on reads on a single chromosome
    total_reads += output_determine_bins[0]; // increase the number of total reads
    Bin[] bins = output_determine_bins[1]; // store the bins
    number_of_bins += bins.length; // increase number of bins
    
    foreach( bin ; bins ) { // loop through the bins
      binsFile.writeln(bin.chr,"\t",bin.start,"\t",bin.end,"\t",bin.size,"\t",bin.fwd,"\t",bin.rev,"\t",bin.cat); // report a single bin
      bin_size_distribution ~= to!double(bin.size); // add the bin size to the bin size distribution
      if (bin.cat == "FR") { FR_bins++; } // increase the number of FR bins
      total_bins++; // increase the total number of bins
    }
    Region[] regions = determine_regions(bins); // determine regions based on bins
    foreach( region ; regions ) { // loop through the regions
      regionsFile.writeln(region.chr,"\t",region.start,"\t",region.end,"\t",region.size,"\t",region.cat); // report a single region
    }    
    Point[] switch_points = determine_switch_points(regions); // determine switch point based on regions
    foreach( point ; switch_points ) { // loop through the switchpoints
      pointsFile.writeln(point.chr,"\t",point.pos); // report a single switch point
    }
  }
  auto FR_bins_ratio = to!double(FR_bins)/to!double(total_bins); // calculcate the ratio of FR bins 
  auto bin_size_median = median(bin_size_distribution); // calucate the median bin size
  auto t = tuple(to!double(total_reads),FR_bins_ratio, bin_size_median,number_of_bins); // combine variables to a single output variable
  return(t); // return output variable
}

// function to determine bins
// A bin is a region with x number of reads in row
Tuple!(double, Bin[]) determine_bins(BamReader bam, ReferenceSequenceInfo chromosome ) {
  // initialize variables
  BamReadBlock[] binreads;
  ulong forward_reads = 0;
  ulong reverse_reads = 0;
  double total_reads = 0;
  Bin[] bins;

  auto reads = bam[chromosome.name]; // get all reads on the single chromsome 
  foreach( read ; reads ) { // loop through the reads
    if (read.mapping_quality >= MQ && !read.is_duplicate) { // skip reads with low mapping_quality and duplicate reads
      if (binreads.length < reads_per_bin) { // if the number of reads in a single bin is smaller then the minimal number of reads in a single bin
	if (read.strand == '+') { forward_reads++; } // increase the number of reads on the forward strand inthe bin
	else if (read.strand == '-') { reverse_reads++; } // increase the number of reads on the reverse strand in the bin
	binreads ~= read; // add a single read to the bin
      } else { // if the number of reads in the bin is equal to the minimal number of reads in a single bin
	string cat = get_bin_category(binreads, forward_reads, reverse_reads); // determine category of the bin (F, R or FR)
	Bin bin = new Bin(chromosome.name, binreads[0].position, binreads[$-1].position, forward_reads, reverse_reads, cat); // create a new bin object
	bins ~= bin; // add the bin to the list of bins
	// if the first read in the bin is on the forward strand and the new read in on the reverse strand
	// decrease the number reads on the forward strand and increase the number of reads on the reverse strand
	if (binreads[0].strand == '+' && read.strand == '-') {
	  forward_reads--;
	  reverse_reads++;
	// if the first read in the bin is on the reverse strand and the new read in on the forward strand
	// decrease the number reads on the reverse strand and increase the number of reads on the forward strand	
	} else if (binreads[0].strand == '-' && read.strand == '+') {
	  reverse_reads--;
	  forward_reads++;
	}
	binreads = binreads[1..$]; // remove first read in the bin from the bin	
	binreads ~= read; // add the new read to the bin	
      }
      total_reads++; // increase the total number of reads
    }
  }
  // process the last bin of the chromsome
  if (binreads.length > 0) {
    string cat = get_bin_category(binreads, forward_reads, reverse_reads); // determine category of the bin (F, R or FR)
    Bin bin = new Bin(chromosome.name, binreads[0].position, binreads[$-1].position, forward_reads, reverse_reads, cat); // create a new bin object
    bins ~= bin; // add the bin to the list of bins
  }
  auto t = tuple(total_reads,bins); // combine variables to a single output variable
  return(t); // return output variable
}

// function to determine category
// A category is based on the number of reads on the forward strand
string get_bin_category(BamReadBlock[] reads, ulong fwd_reads, ulong rev_reads) {
  string category = "FR"; // initialize category variable  
  double fwd_ratio = to!double(fwd_reads)/(to!double(rev_reads)+to!double(fwd_reads)); // calculate 'forward' reads ratio
  if (fwd_ratio >= 0.99) { category = "F"; } // category = F(orward) if 99% of the reads is on the forward strand
  else if (fwd_ratio <= 0.01) { category = "R"; } // category = R(everse) if 1% of the reads in on the forward strand
  return(category); // return category variable
}

// function to determing regions
// A region is an x number of bins in a row with the same category
Region[] determine_regions(Bin[] bins) {
  // initialize variables
  string prev_cat = "X";
  Region[] regions_out;
  Bin[][string] regions;
  
  foreach( bin ; bins ) { // loop through each bin
    if (prev_cat != "X") { // if it's not the first  bin
      if (bin.cat == prev_cat) { // if category is equal to the category of the previous bin
	foreach( region_cat ; regions.keys ) { // loop through the different categories in the region
	  // if the number of bins with the same category is bigger or equal to the minimal number of reads in a bin
	  // and the category is different then the current category
	  // and the number of bins with the current category in the region is bigger or equal to the minimal number of reads in a bin
	  if (regions[region_cat].length >= reads_per_bin && region_cat != bin.cat && regions[bin.cat].length >= reads_per_bin) {
	    Region region = new Region(regions[region_cat][0].chr, regions[region_cat][0].start, regions[region_cat][$-1].end, region_cat); // create a new region object
	    regions_out ~= region; // add region to a list of regions
	    regions[region_cat] = null; // empty the region with the specific category
	    // if number of bins in region is lower then the minimal number of reads in a bin and the category is different then the current category
	  } else if (regions[region_cat].length < reads_per_bin && region_cat != bin.cat) { 
	    regions[region_cat] = null; // empty the region with the specific category
	  }
	}
	regions[bin.cat] ~= bin; // add bin to the region with the specific category
      } else { // if category is different then the previous category
	// if current category exists in the region and the number of bins is bigger or equal to the minimal number of reads in a bin
	if (bin.cat in regions && regions[bin.cat].length >= reads_per_bin) { 
	  foreach( region_cat ; regions.keys ) { // loop through all categories in the region
	    if (region_cat != bin.cat) {  regions[region_cat] = null; } // if category is different then the current category, discard category from region
	  }
	}
	regions[bin.cat] ~= bin; // add bin to the region with the current category
      }
    } else { regions[bin.cat] ~= bin; } // add the first bin to the region with the current category
    prev_cat = bin.cat.dup; // store current category to the previous category
  }
  foreach( region_cat ; regions.keys ) { // loop trough the categories in the last region
    if (regions[region_cat].length > 1) { // if the number of bins is bigger the 1
      Region region = new Region(regions[region_cat][0].chr, regions[region_cat][0].start, regions[region_cat][$-1].end, region_cat); // create new region object
      regions_out ~= region; // add region to a list of region
    }
  }
  return(regions_out); // return the list of regions
}

// function to determine switch points
// A switch point is the start position of a region
Point[] determine_switch_points(Region[] regions) {
  Point[] switch_points; // initialize list of switch points
  for( int i = 1; i < regions.length; i++) {
//   foreach( region ; regions ) { // loop through the regions
    Point point = new Point(regions[i].chr, regions[i].start); // create a switch point object
    switch_points ~= point; // add the switch point to the list of switch points
  }
  return(switch_points); // return the list of switch point
}

// function for the QC analysis
double[] qc_analysis(DirEntry[] bamFiles, double[] total_reads_distribution, double[] FR_bins_ratio_distribution, double[] bin_size_median_distribution, double[] number_of_bins_distribution, string folder) {
  // initialize variables
  double[] bin_size_median_distribution2;
  ulong passed_cells = 0;
    
  // calculate means and standard deviations
  double total_reads_mean = mean(total_reads_distribution);
  double total_reads_sd = stdev(total_reads_distribution);
  double FR_bins_ratio_mean = mean(FR_bins_ratio_distribution);
  double FR_bins_ratio_sd = stdev(FR_bins_ratio_distribution);    
  double number_of_bins_mean = mean(number_of_bins_distribution);
  double number_of_bins_sd = stdev(number_of_bins_distribution);
  
  string outName = get_output_name(folder); // get the output name
  File statsFile = File(outName~"_stats.txt", "w"); // open the output file for the statistics
  
  // write the header to the statistics output file
  statsFile.writeln("##TOTAL READS MEAN: ",total_reads_mean);
  statsFile.writeln("##TOTAL READS SD: ",total_reads_sd);
  statsFile.writeln("##FR BIN RATIO MEAN: ",FR_bins_ratio_mean);
  statsFile.writeln("##FR BIN RATIO SD: ",FR_bins_ratio_sd);
  statsFile.writeln("##NUMBER OF BINS MEAN: ", number_of_bins_mean);
  statsFile.writeln("##NUMBER OF BINS SD: ",number_of_bins_sd);
  statsFile.writeln("#FILE\tTOTAL_READS\tFR_BIN_RATIO\tMEDIAN_BIN_SIZE\tREADS_PER_BINS\tNUMBER_OF_BINS\tFILTER\tREASON");
  
  for( int i; i < bamFiles.length; i++ ) { // loop through all the *.bam files
    // if the total number of reads is less then the mean-sd, this bam file failed the QC
    if (total_reads_distribution[i] < (total_reads_mean-total_reads_sd)) {
      qc_fail[folder][bamFiles[i]~".points"] = 1; // store qc failed cells
      statsFile.writeln(bamFiles[i],"\t",total_reads_distribution[i],"\t",FR_bins_ratio_distribution[i],"\t",bin_size_median_distribution[i],"\t",to!int(reads_per_bin),"\t",number_of_bins_distribution[i],"\t","FAILED","\t","TOTAL READS");
    // if the ratio of FR bins is bigger then the mean+sd, this bam file failed the QC
    } else if (FR_bins_ratio_distribution[i] > (FR_bins_ratio_mean+FR_bins_ratio_sd)) {
      qc_fail[folder][bamFiles[i]~".points"] = 1; // store qc failed cells
      statsFile.writeln(bamFiles[i],"\t",total_reads_distribution[i],"\t",FR_bins_ratio_distribution[i],"\t",bin_size_median_distribution[i],"\t",to!int(reads_per_bin),"\t",number_of_bins_distribution[i],"\t","FAILED","\t","FR BIN RATIO");
    // if the total number of bins is less then the mean+sd, this bam file failed the QC    
    } else if (number_of_bins_distribution[i] < (number_of_bins_mean-number_of_bins_sd)) {
      qc_fail[folder][bamFiles[i]~".points"] = 1; // store qc failed cells
      statsFile.writeln(bamFiles[i],"\t",total_reads_distribution[i],"\t",FR_bins_ratio_distribution[i],"\t",bin_size_median_distribution[i],"\t",to!int(reads_per_bin),"\t",number_of_bins_distribution[i],"\t","FAILED","\t","NUMBER OF BINS");
    // report the passed QC bam files
    } else {
      statsFile.writeln(bamFiles[i],"\t",total_reads_distribution[i],"\t",FR_bins_ratio_distribution[i],"\t",bin_size_median_distribution[i],"\t",to!int(reads_per_bin),"\t",number_of_bins_distribution[i],"\t","PASS","\t","");
      bin_size_median_distribution2 ~= bin_size_median_distribution[i]; // add median distribution of a passed QC bam file
      passed_cells++; // increase the number of passed QC files
    }
  }
  // calculate the minimal number of cells, based on only passed QC files
  if (outName == patientName) { patient_count = to!double(passed_cells)/100*to!double(patient_count_perc); } 
  else if (outName == fatherName) { father_count = to!double(passed_cells)/100*to!double(patient_count_perc); } 
  else if (outName == motherName) { mother_count = to!double(passed_cells)/100*to!double(patient_count_perc); }
  return(bin_size_median_distribution2); // return meadian distribution of only passed QC files
}

// function to merge switch points
void merge_switch_points(string patientFolder, string fatherFolder, string motherFolder, double merge_distance) {
  
  ulong[DirEntry][string][ulong][string] switch_points; // initialize switch points variable  
  
  // initialize output file name
  string outputName = "switch_points_p"~patientName~"_f"~fatherName~"_m"~motherName~"_r"~to!string(reads_per_bin)~"_d"~to!string(distance_factor)~"_x"~to!string(patient_percentage)~"_c"~to!string(patient_count_perc)~"_q"~to!string(MQ)~".txt";
  File output = File(outputName,"w"); // open output file
  
  // write header to the output file
  output.writeln("##PATIENT: ",patientName);
  output.writeln("##FATHER: ",fatherName);
  output.writeln("##MOTHER: ",motherName);
  output.writeln("##READS: ",reads_per_bin);
  output.writeln("##DISTANCE FACTOR: ",distance_factor);
  output.writeln("##PERCENTAGE PATIENT CELLS: ", patient_percentage);
  output.writeln("##PATIENT COUNT PERCENTAGE: ",patient_count_perc);
  output.writeln("##MERGE_DISTANCE: ",merge_distance);  
  output.writeln("##PATIENT COUNT: ",to!ulong(patient_count));  
  output.writeln("#CHR\tSTART\tEND\t",patientName,"\t",fatherName,"\t",motherName,"\tPATIENT_PERC\tFATHER_PERC\tMOTHER_PERC\tFILTER\tREASON\tPATIENT_UNIQ\tFATHER_UNIQ\tMOTHER_UNIQ");
    
  foreach( folder ; [patientFolder, fatherFolder, motherFolder] ) { // loop trough each folder of all samples
    string name = get_output_name(folder); // get output name of the sample
    auto pointsFiles = array(filter!`endsWith(a.name,".points")`(dirEntries(folder,SpanMode.depth))); // get all *.points files
    foreach( pointsFile ; pointsFiles ) { // loop through all switch points files      
      if (pointsFile in qc_fail[folder]) { continue; } // skip qc failed cells
      File file = File(pointsFile,"r"); // open switch point file
      while(!file.eof()) { // loop through the switch point file
	string[] line = chomp(file.readln()).split("\t"); // split each line on tabs
	if (line.length > 1 && line[0] != "#CHR") { switch_points[line[0]][to!ulong(line[1])][name][pointsFile]++; } // store the switch point in a hash of switch points
      }
    }
  }
  
  foreach( chr ; switch_points.keys.sort ) { // loop through all switch points in all files of all samples
    // initialize variables
    ulong prev_pos = -1;
    ulong[] switch_point_region;
    ulong[string][string] samples;
    
    foreach( pos ; switch_points[chr].keys.sort ) { // loop through each switch point position on a specific chromosome
      foreach( sample ; switch_points[chr][pos].keys ) { // loop through each sample containing the switch point position
	foreach( file ; switch_points[chr][pos][sample].keys ) { // loop through each file of the sample containing the switch point position
	  if ((pos-prev_pos) < merge_distance || prev_pos == -1) { // if the distance between two switch points is less then 'variable', except the first switch point position
	    switch_point_region ~= pos; // store position in a switch point region
	    samples[sample][file]++; // store the sample and file in a switch point region
	  } else if (prev_pos != -1) { // if the distance bewtween two switch_points is bigger or equal to 'variable', except the first switch point position
	    // initialize sample specific variables
	    ulong patient_switch_points = 0;
	    ulong father_switch_points = 0;
	    ulong mother_switch_points = 0;
	    int patient_uniq = 0;
	    int father_uniq = 0;
	    int mother_uniq = 0;	    
	    
	    // increase the number of switch point in the switch point region found in each sample
	    if (patientName in samples) { patient_switch_points = samples[patientName].keys.length; }
	    if (fatherName in samples) { father_switch_points = samples[fatherName].keys.length; }
	    if (motherName in samples) { mother_switch_points = samples[motherName].keys.length; }
	    
	    // calculate the percentage of switch points in the switch point region found in each sample
	    double perc_patient = patient_switch_points/to!double(patient_switch_points+father_switch_points+mother_switch_points)*100;
	    double perc_father = father_switch_points/to!double(patient_switch_points+father_switch_points+mother_switch_points)*100;
	    double perc_mother = mother_switch_points/to!double(patient_switch_points+father_switch_points+mother_switch_points)*100;
	    
	    // determine sample specific switch point regions
	    if (perc_patient >= patient_percentage && patient_switch_points >= patient_count) { patient_uniq = 1; }
	    else if (perc_father >= patient_percentage && father_switch_points >= father_count) { father_uniq = 1; }
	    else if (perc_mother >= patient_percentage && mother_switch_points >= mother_count) { mother_uniq = 1; }
	    
	    // if the percentage of switch points in the switch point region found in the patient is less then 'variable', report as failed
	    if (perc_patient < patient_percentage) {
	      output.writeln(chr,"\t",switch_point_region[0],"\t",switch_point_region[$-1],"\t",patient_switch_points,"\t",father_switch_points,"\t",mother_switch_points,"\t",perc_patient,"\t",perc_father,"\t",perc_mother,"\t","FAILED","\t","PATIENT PERCENTAGE","\t",patient_uniq,"\t",father_uniq,"\t",mother_uniq);
	    // if the number of switch points in the switch point region found in the patient is less then 'variable', report as failed	      
	    } else if (patient_switch_points < patient_count) {
	      output.writeln(chr,"\t",switch_point_region[0],"\t",switch_point_region[$-1],"\t",patient_switch_points,"\t",father_switch_points,"\t",mother_switch_points,"\t",perc_patient,"\t",perc_father,"\t",perc_mother,"\t","FAILED","\t","PATIENT COUNT","\t",patient_uniq,"\t",father_uniq,"\t",mother_uniq);
	    // report patient denovo switch points
	    } else {
	      output.writeln(chr,"\t",switch_point_region[0],"\t",switch_point_region[$-1],"\t",patient_switch_points,"\t",father_switch_points,"\t",mother_switch_points, "\t",perc_patient,"\t",perc_father,"\t",perc_mother,"\t","PASS","\t","-","\t",patient_uniq,"\t",father_uniq,"\t",mother_uniq);
	    }
	    switch_point_region = [pos]; // empty switch point region and store first new switch point position
	    samples = null; // empty samples list
	    samples[sample][file]++; // add sample and file to the list
	  }
	}
      }
      prev_pos = pos; // store switch point position to the previous switch point position
    }
    // last switch point region
    // initialize sample specific variables
    ulong patient_switch_points = 0;
    ulong father_switch_points = 0;
    ulong mother_switch_points = 0;
    int patient_uniq = 0;
    int father_uniq = 0;
    int mother_uniq = 0;	    
    
    // increase the number of switch point in the switch point region found in each sample
    if (patientName in samples) { patient_switch_points = samples[patientName].keys.length; }
    if (fatherName in samples) { father_switch_points = samples[fatherName].keys.length; }
    if (motherName in samples) { mother_switch_points = samples[motherName].keys.length; }
    
    // calculate the percentage of switch points in the switch point region found in each sample
    double perc_patient = patient_switch_points/to!double(patient_switch_points+father_switch_points+mother_switch_points)*100;
    double perc_father = father_switch_points/to!double(patient_switch_points+father_switch_points+mother_switch_points)*100;
    double perc_mother = mother_switch_points/to!double(patient_switch_points+father_switch_points+mother_switch_points)*100;
    
    // determine sample specific switch point regions
    if (perc_patient >= patient_percentage && patient_switch_points >= patient_count) { patient_uniq = 1; }
    else if (perc_father >= patient_percentage && father_switch_points >= father_count) { father_uniq = 1; }
    else if (perc_mother >= patient_percentage && mother_switch_points >= mother_count) { mother_uniq = 1; }
    
    // if the percentage of switch points in the switch point region found in the patient is less then 'variable', report as failed
    if (perc_patient < patient_percentage) {
      output.writeln(chr,"\t",switch_point_region[0],"\t",switch_point_region[$-1],"\t",patient_switch_points,"\t",father_switch_points,"\t",mother_switch_points,"\t",perc_patient,"\t",perc_father,"\t",perc_mother,"\t","FAILED","\t","PATIENT PERCENTAGE","\t",patient_uniq,"\t",father_uniq,"\t",mother_uniq);
    // if the number of switch points in the switch point region found in the patient is less then 'variable', report as failed	      
    } else if (patient_switch_points < patient_count) {
      output.writeln(chr,"\t",switch_point_region[0],"\t",switch_point_region[$-1],"\t",patient_switch_points,"\t",father_switch_points,"\t",mother_switch_points,"\t",perc_patient,"\t",perc_father,"\t",perc_mother,"\t","FAILED","\t","PATIENT COUNT","\t",patient_uniq,"\t",father_uniq,"\t",mother_uniq);
    // report patient denovo switch points
    } else {
      output.writeln(chr,"\t",switch_point_region[0],"\t",switch_point_region[$-1],"\t",patient_switch_points,"\t",father_switch_points,"\t",mother_switch_points, "\t",perc_patient,"\t",perc_father,"\t",perc_mother,"\t","PASS","\t","-","\t",patient_uniq,"\t",father_uniq,"\t",mother_uniq);
    }
  }  
}

// function to get the output name
string get_output_name( string folder ) {
  if (folder.endsWith("/")) { folder = replace(folder,regex(r"\/$"),""); } // parse last '/' character of the input folder name
  string[] path = folder.split("/"); // split input folder name on '/'
  string outName = path[$-1]; // store last part of the input folder name
  return outName; // return output name
}

// function to calculate the mean of an array
double mean(double[] array) {
  double sum = 0.0; // initialize sum variable
  foreach( item ; array ) { // loop throuh each item of the array
    sum += item; // increase the sum
  }
  double mean = sum/to!double(array.length); // calculate the mean
  return(mean); // return the mean
}

// function to calculate the median of an array
double median(double[] array) {
  double median = 0; // initialize median variable  
  array.sort; // sort array
  double pos = to!double(array.length/2); // get middle position of the array
  if ((array.length % 2) != 0) { median = array[to!ulong(pos+0.5)]; } // if the length of the array is odd, median is the middle value
  else { median = (array[to!ulong(pos)]+array[to!ulong(pos-1)])/2; } // if the length of the array is even, median is the mean of the two middle values
  return(median); // return median
}

// function to calculate the standard deviation of an array
double stdev(double[] array) {
  double sqtotal = 0.0; // initialize the square total variable
  double mean = mean(array); // calculate the man of the array
  foreach( item ; array ) { // loop through each item of the array
    sqtotal += (mean-to!double(item))^^2; // calculate the square root distance of the variable from the mean
  }
  auto std = (sqtotal/array.length)^^0.5; // calculate the standard deviation
  return(std); // return standard deviation
}

// class to create bin object
class Bin {
  string chr;
  ulong start;
  ulong end;
  ulong size;
  ulong fwd;
  ulong rev;
  string cat;  
  this(string chr, ulong start, ulong end, ulong fwd, ulong rev, string cat) {
    this.chr = chr;
    this.start = start;
    this.end = end;
    this.size = abs(end-start);
    this.fwd = fwd;
    this.rev = rev;
    this.cat = cat;
  }
}

// class to create region object
class Region {
  string chr;
  ulong start;
  ulong end;
  ulong size;
  string cat;  
  this(string chr, ulong start, ulong end, string cat) {
    this.chr = chr;
    this.start = start;
    this.end = end;
    this.size = abs(end-start);
    this.cat = cat;
  }
}

// class to create switch point object
class Point {
  string chr;
  ulong pos;
  this(string chr, ulong pos) {
    this.chr = chr;
    this.pos = pos;
  }
}
