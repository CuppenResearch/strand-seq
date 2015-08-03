 
// Sliding dynamic window binning.
// Remove duplicates
// 1% category


import std.getopt, std.math, std.stdio, std.regex, std.algorithm, std.range, std.random, std.file, std.string, std.conv;
import bio.bam.reader, bio.bam.pileup;

bool help;
string patientFolder;
string fatherFolder;
string motherFolder;
int MQ = 29;
ulong reads_per_bin = 30;
double distance_factor = 2;
double patient_percentage = 65;
double patient_count_perc = 25;

double patient_count;
double father_count;
double mother_count;
string patientName;
string fatherName;
string motherName;

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
    stderr.writeln("                    number of reads per bin (DEFAULT: 30)");
    stderr.writeln("         -d, --distance");
    stderr.writeln("                    distance factor (DEFAULT: 2)");
    stderr.writeln("         -x, --percentage");
    stderr.writeln("                    minimal percentage of patient cells of all cells supporting the switchpoint(DEFAULT: 65%)");
    stderr.writeln("         -c, --count");
    stderr.writeln("                    minimal percentage of 'passed' cells in the patient supporting the switchpoint (DEFAULT: 25%)");
    stderr.writeln("         -q, --MQ");
    stderr.writeln("                    skip reads with mapping quality le value (DEFAULT: 10)");
}

int main(string[] args) {
    getopt(args, 
    std.getopt.config.caseSensitive,
    std.getopt.config.passThrough,
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
  if (help || patientFolder == null || fatherFolder == null || motherFolder == null) {
    printUsage();
    return 0;
  }
  patientName = get_output_name(patientFolder);
  fatherName = get_output_name(fatherFolder);
  motherName = get_output_name(motherFolder);
  
  double[] bin_size_median_distribution_patient = bam_looping(patientFolder);
  double[] bin_size_median_distribution_father = bam_looping(fatherFolder);
  double[] bin_size_median_distribution_mother = bam_looping(motherFolder);
  
  double[] bin_size_median_distribution_trio = bin_size_median_distribution_patient ~ bin_size_median_distribution_father ~ bin_size_median_distribution_mother;
  double bin_size_median_trio_mean = mean(bin_size_median_distribution_trio);
  double merge_distance = bin_size_median_trio_mean/distance_factor;
  
  merge_switch_points(patientFolder, fatherFolder, motherFolder, merge_distance);  
  
  return 1;
}

double[] bam_looping( string folder ) {
  auto bamFiles = array(filter!`endsWith(a.name,".bam")`(dirEntries(folder,SpanMode.depth)));
  double[] total_reads_distribution;
  double[] FR_bins_ratio_distribution;
  double[] bin_size_median_distribution;
  double[] number_of_bins_distribution;
  foreach( bamFile ; bamFiles ) {
    auto output = process_bam(bamFile);
    total_reads_distribution ~= output[0];
    FR_bins_ratio_distribution ~= output[1];
    bin_size_median_distribution ~= output[2];
    number_of_bins_distribution ~= output[3];
  }  
  bin_size_median_distribution = qc_analysis(bamFiles, total_reads_distribution, FR_bins_ratio_distribution, bin_size_median_distribution, number_of_bins_distribution, folder);  
 
  return(bin_size_median_distribution);
}

Tuple!(double, double, double, double) process_bam(string file) {
  stdout.writeln("Busy with "~file~"...");
  BamReader bam = new BamReader(file);
  double total_reads = 0;
  double number_of_bins = 0;
  double[] bin_size_distribution;
  ulong FR_bins;
  ulong total_bins;    
  File binsFile = File(file~".bins","w");
  File regionsFile = File(file~".regions","w");
  File pointsFile = File(file~".points","w");
  
  binsFile.writeln("#CHR\tSTART\tEND\tSIZE\tF\tR\tCAT");
  regionsFile.writeln("#CHR\tSTART\tEND\tSIZE\tCAT");
  pointsFile.writeln("#CHR\tPOS");
  
  foreach( chromosome ;  bam.reference_sequences() ) {
    if (startsWith(chromosome.name, "GL")) {
      continue;
    }
    auto output_determine_bins = determine_bins(bam, chromosome);
    total_reads += output_determine_bins[0];
    Bin[] bins = output_determine_bins[1];
    number_of_bins += bins.length;
    foreach( bin ; bins ) {
      binsFile.writeln(bin.chr,"\t",bin.start,"\t",bin.end,"\t",bin.size,"\t",bin.fwd,"\t",bin.rev,"\t",bin.cat);
      bin_size_distribution ~= to!double(bin.size);
      if (bin.cat == "FR") {
	FR_bins++;
      }
      total_bins++;
    }
    Region[] regions = determine_regions(bins);
    foreach( region ; regions ) {
      regionsFile.writeln(region.chr,"\t",region.start,"\t",region.end,"\t",region.size,"\t",region.cat);
    }    
    Point[] switch_points = determine_switch_points(regions);
    foreach( point ; switch_points ) {
      pointsFile.writeln(point.chr,"\t",point.pos);
    }
  }
  auto FR_bins_ratio = to!double(FR_bins)/to!double(total_bins);
  auto bin_size_median = median(bin_size_distribution);
  auto t = tuple(to!double(total_reads),FR_bins_ratio, bin_size_median,number_of_bins);
  return(t);
}

Tuple!(double, Bin[]) determine_bins(BamReader bam, ReferenceSequenceInfo chromosome ) {
  auto reads = bam[chromosome.name];
  BamReadBlock[] binreads;
  ulong forward_reads = 0;
  ulong reverse_reads = 0;
  double total_reads = 0;
  Bin[] bins;

  foreach( read ; reads ) {
    if (read.mapping_quality > MQ && !read.is_duplicate) {
      if (binreads.length < reads_per_bin) {
	if (read.strand == '+') {
	  forward_reads++;
	} else if (read.strand == '-') {
	  reverse_reads++;
	}
	binreads ~= read;
      } else {
	string cat = get_bin_category(binreads, forward_reads, reverse_reads);
	Bin bin = new Bin(chromosome.name, binreads[0].position, binreads[$-1].position, forward_reads, reverse_reads, cat);
	bins ~= bin;	  
	if (binreads[0].strand == '+' && read.strand == '-') {
	  forward_reads--;
	  reverse_reads++;
	} else if (binreads[0].strand == '-' && read.strand == '+') {
	  reverse_reads--;
	  forward_reads++;
	}
	binreads = binreads[1..$];
	binreads ~= read;
      }
      total_reads++;	
    }
  }
  if (binreads.length > 0) {
    string cat = get_bin_category(binreads, forward_reads, reverse_reads);
    Bin bin = new Bin(chromosome.name, binreads[0].position, binreads[$-1].position, forward_reads, reverse_reads, cat);
    bins ~= bin;
  }
  auto t = tuple(total_reads,bins);
  return(t);
}

string get_bin_category(BamReadBlock[] reads, ulong fwd_reads, ulong rev_reads) {
  int start = reads[0].position;
  int end = reads[$-1].position;
  double fwd_ratio = to!double(fwd_reads)/(to!double(rev_reads)+to!double(fwd_reads));
  string category = "FR";
  if (fwd_ratio >= 0.9) {
    category = "F";
  } else if (fwd_ratio <= 0.1) {
    category = "R";
  }
  return(category);
}

Region[] determine_regions(Bin[] bins) {
  string prev_cat = "X";
  Region[] regions_out;
  Bin[][string] regions;
  foreach( bin ; bins ) {
    if (prev_cat != "X") {
      if (bin.cat == prev_cat) {
	foreach( region_cat ; regions.keys ) {
	  if (regions[region_cat].length >= reads_per_bin && region_cat != bin.cat && regions[bin.cat].length >= reads_per_bin) {
	    Region region = new Region(regions[region_cat][0].chr, regions[region_cat][0].start, regions[region_cat][$-1].end, region_cat);
	    regions_out ~= region;
	    regions[region_cat] = null;
	  } else if (regions[region_cat].length < reads_per_bin && region_cat != bin.cat) {
	    regions[region_cat] = null;
	  }
	}
	regions[bin.cat] ~= bin;
      } else {
	if (bin.cat in regions && regions[bin.cat].length >= reads_per_bin) {
	  foreach( region_cat ; regions.keys ) {
	    if (region_cat != bin.cat) {
	      regions[region_cat] = null;
	    }
	  }
	}
	regions[bin.cat] ~= bin;
      }
    } else {
      regions[bin.cat] ~= bin;
    }	  
    prev_cat = bin.cat.dup;
  }
  foreach( region_cat ; regions.keys ) {
    if (regions[region_cat].length > 1) {
      Region region = new Region(regions[region_cat][0].chr, regions[region_cat][0].start, regions[region_cat][$-1].end, region_cat);
      regions_out ~= region;
    }
  }
  return(regions_out);
}

Point[] determine_switch_points(Region[] regions) {
  Point[] switch_points;
  for( int i = 1; i < regions.length; i++ ) {
    Point point = new Point(regions[i].chr, regions[i].start);
    switch_points ~= point;
  }
  return(switch_points);
}

double[] qc_analysis(DirEntry[] bamFiles, double[] total_reads_distribution, double[] FR_bins_ratio_distribution, double[] bin_size_median_distribution, double[] number_of_bins_distribution, string folder) {
  string[] failed;  
  double total_reads_mean = mean(total_reads_distribution);
  double total_reads_sd = stdev(total_reads_distribution);
  double FR_bins_ratio_mean = mean(FR_bins_ratio_distribution);
  double FR_bins_ratio_sd = stdev(FR_bins_ratio_distribution);    
  double number_of_bins_mean = mean(number_of_bins_distribution);
  double number_of_bins_sd = stdev(number_of_bins_distribution);
  double[] bin_size_median_distribution2;
  string outName = get_output_name(folder);  
  ulong passed_cells = 0;
  File statsFile = File(outName~"_stats.txt", "w");
  statsFile.writeln("##TOTAL READS MEAN: ",total_reads_mean);
  statsFile.writeln("##TOTAL READS SD: ",total_reads_sd);
  statsFile.writeln("##FR BIN RATIO MEAN: ",FR_bins_ratio_mean);
  statsFile.writeln("##FR BIN RATIO SD: ",FR_bins_ratio_sd);
  statsFile.writeln("##NUMBER OF BINS MEAN: ", number_of_bins_mean);
  statsFile.writeln("##NUMBER OF BINS SD: ",number_of_bins_sd);
  statsFile.writeln("#FILE\tTOTAL_READS\tFR_BIN_RATIO\tMEDIAN_BIN_SIZE\tREADS_PER_BINS\tNUMBER_OF_BINS\tFILTER\tREASON");
  for( int i; i < bamFiles.length; i++ ) {
//     auto reads_per_bin = total_reads_distribution[i]/100*reads_percentage;
    if (total_reads_distribution[i] < (total_reads_mean-total_reads_sd)) {
      statsFile.writeln(bamFiles[i],"\t",total_reads_distribution[i],"\t",FR_bins_ratio_distribution[i],"\t",bin_size_median_distribution[i],"\t",to!int(reads_per_bin),"\t",number_of_bins_distribution[i],"\t","FAILED","\t","TOTAL READS");
    } else if (FR_bins_ratio_distribution[i] > (FR_bins_ratio_mean+FR_bins_ratio_sd)) {
      statsFile.writeln(bamFiles[i],"\t",total_reads_distribution[i],"\t",FR_bins_ratio_distribution[i],"\t",bin_size_median_distribution[i],"\t",to!int(reads_per_bin),"\t",number_of_bins_distribution[i],"\t","FAILED","\t","FR BIN RATIO");
    } else if (reads_per_bin < 2) {
      statsFile.writeln(bamFiles[i],"\t",total_reads_distribution[i],"\t",FR_bins_ratio_distribution[i],"\t",bin_size_median_distribution[i],"\t",to!int(reads_per_bin),"\t",number_of_bins_distribution[i],"\t","FAILED","\t","READS PER BIN");
    } else if (number_of_bins_distribution[i] < (number_of_bins_mean-number_of_bins_sd)) {
      statsFile.writeln(bamFiles[i],"\t",total_reads_distribution[i],"\t",FR_bins_ratio_distribution[i],"\t",bin_size_median_distribution[i],"\t",to!int(reads_per_bin),"\t",number_of_bins_distribution[i],"\t","FAILED","\t","NUMBER OF BINS");
    } else {
      statsFile.writeln(bamFiles[i],"\t",total_reads_distribution[i],"\t",FR_bins_ratio_distribution[i],"\t",bin_size_median_distribution[i],"\t",to!int(reads_per_bin),"\t",number_of_bins_distribution[i],"\t","PASS","\t","");
      bin_size_median_distribution2 ~= bin_size_median_distribution[i];
      passed_cells++;
    }
  }
  if (outName == patientName) {
    patient_count = to!double(passed_cells)/100*to!double(patient_count_perc);
  } else if (outName == fatherName) {
    father_count = to!double(passed_cells)/100*to!double(patient_count_perc);
  } else if (outName == motherName) {
    mother_count = to!double(passed_cells)/100*to!double(patient_count_perc);
  }
  return(bin_size_median_distribution2);
}




void merge_switch_points(string patientFolder, string fatherFolder, string motherFolder, double merge_distance) {
  
  string outputName = "switch_points_p"~patientName~"_f"~fatherName~"_m"~motherName~"_r"~to!string(reads_per_bin)~"_d"~to!string(distance_factor)~"_x"~to!string(patient_percentage)~"_c"~to!string(patient_count_perc)~"_q"~to!string(MQ)~".txt";
    
  File output = File(outputName,"w");
  
  output.writeln("##PATIENT: ",patientName);
  output.writeln("##FATHER: ",fatherName);
  output.writeln("##MOTHER: ",motherName);
  output.writeln("##READS: ",reads_per_bin);
  output.writeln("##DISTANCE FACTOR: ",distance_factor);
  output.writeln("##PERCENTAGE PATIENT CELLS: ", patient_percentage);
  output.writeln("##PATIENT COUNT PERCENTAGE: ",patient_count_perc);
  output.writeln("##MERGE_DISTANCE: ",merge_distance);  
  output.writeln("##PATIENT COUNT: ",to!ulong(patient_count));
  
  output.writeln("#CHR\tSTART\tEND\t","TUNED","\t",patientName,"\t",fatherName,"\t",motherName,"\tPATIENT_PERC\tFATHER_PERC\tMOTHER_PERC\tFILTER\tREASON\tPATIENT_UNIQ\tFATHER_UNIQ\tMOTHER_UNIQ");
  
  ulong[DirEntry][string][ulong][string] switch_points;
  
  foreach( folder ; [patientFolder, fatherFolder, motherFolder] ) {
    string name = get_output_name(folder);
    auto pointsFiles = array(filter!`endsWith(a.name,".points")`(dirEntries(folder,SpanMode.depth)));
    foreach( pointsFile ; pointsFiles ) {
      File file = File(pointsFile,"r");
      while(!file.eof()) {
	string[] line = chomp(file.readln()).split("\t");
	if (line.length > 1 && line[0] != "#CHR") {
	  switch_points[line[0]][to!ulong(line[1])][name][pointsFile]++;
	}
      }
    }
  }
  
  foreach( chr ; switch_points.keys.sort ) {
    ulong prev_pos = -1;
    ulong[] switch_point_region;
    ulong[string][string] samples;
    foreach( pos ; switch_points[chr].keys.sort ) {
      foreach( sample ; switch_points[chr][pos].keys ) {
	foreach( file ; switch_points[chr][pos][sample].keys ) {
	  if ((pos-prev_pos) < merge_distance || prev_pos == -1) {
	    switch_point_region ~= pos;
	    samples[sample][file]++;
	  } else if (prev_pos != -1) {
	    ulong patient_switch_points = 0;
	    ulong father_switch_points = 0;
	    ulong mother_switch_points = 0;
	    if (patientName in samples) {
	      patient_switch_points = samples[patientName].keys.length;
	    }
	    if (fatherName in samples) {
	      father_switch_points = samples[fatherName].keys.length;
	    }
	    if (motherName in samples) {
	      mother_switch_points = samples[motherName].keys.length;
	    }
	    double perc_patient = patient_switch_points/to!double(patient_switch_points+father_switch_points+mother_switch_points)*100;
	    double perc_father = father_switch_points/to!double(patient_switch_points+father_switch_points+mother_switch_points)*100;
	    double perc_mother = mother_switch_points/to!double(patient_switch_points+father_switch_points+mother_switch_points)*100;
	    int patient_uniq = 0;
	    int father_uniq = 0;
	    int mother_uniq = 0;
	    if (perc_patient >= patient_percentage && patient_switch_points >= patient_count) {
	      patient_uniq = 1;
	    } else if (perc_father >= patient_percentage && father_switch_points >= father_count) {
	      father_uniq = 1;
	    } else if (perc_mother >= patient_percentage && mother_switch_points >= mother_count) {
	      mother_uniq = 1;
	    }
	    if (perc_patient < patient_percentage) {
	      output.writeln(chr,"\t",switch_point_region[0],"\t",switch_point_region[$-1],"\t","-","\t",patient_switch_points,"\t",father_switch_points,"\t",mother_switch_points,"\t",perc_patient,"\t",perc_father,"\t",perc_mother,"\t","FAILED","\t","PATIENT PERCENTAGE","\t",patient_uniq,"\t",father_uniq,"\t",mother_uniq);
	    } else if (patient_switch_points < patient_count) {
	      output.writeln(chr,"\t",switch_point_region[0],"\t",switch_point_region[$-1],"\t","-","\t",patient_switch_points,"\t",father_switch_points,"\t",mother_switch_points,"\t",perc_patient,"\t",perc_father,"\t",perc_mother,"\t","FAILED","\t","PATIENT COUNT","\t",patient_uniq,"\t",father_uniq,"\t",mother_uniq);
	    } else {
// 	      double tuned = fine_tuning(chr,switch_point_region[0], switch_point_region[$-1], samples[patientName].keys);
	      double tuned = 0.0;
	      output.writeln(chr,"\t",switch_point_region[0],"\t",switch_point_region[$-1],"\t",to!ulong(tuned),"\t",patient_switch_points,"\t",father_switch_points,"\t",mother_switch_points, "\t",perc_patient,"\t",perc_father,"\t",perc_mother,"\t","PASS","\t","-","\t",patient_uniq,"\t",father_uniq,"\t",mother_uniq);
	    }
	    switch_point_region = [pos];
	    samples = null;
	    samples[sample][file]++;
	  }
	}
      }
      prev_pos = pos;
    }
    ulong patient_switch_points = 0;
    ulong father_switch_points = 0;
    ulong mother_switch_points = 0;
    if (patientName in samples) {
      patient_switch_points = samples[patientName].keys.length;
    }
    if (fatherName in samples) {
      father_switch_points = samples[fatherName].keys.length;
    }
    if (motherName in samples) {
      mother_switch_points = samples[motherName].keys.length;
    }
    double perc_patient = patient_switch_points/to!double(patient_switch_points+father_switch_points+mother_switch_points)*100;
    double perc_father = father_switch_points/to!double(patient_switch_points+father_switch_points+mother_switch_points)*100;
    double perc_mother = mother_switch_points/to!double(patient_switch_points+father_switch_points+mother_switch_points)*100;
    int patient_uniq = 0;
    int father_uniq = 0;
    int mother_uniq = 0;
    if (perc_patient >= patient_percentage && patient_switch_points >= patient_count) {
      patient_uniq = 1;
    } else if (perc_father >= patient_percentage && father_switch_points >= father_count) {
      father_uniq = 1;
    } else if (perc_mother >= patient_percentage && mother_switch_points >= mother_count) {
      mother_uniq = 1;
    }
    if (perc_patient < patient_percentage) {
      output.writeln(chr,"\t",switch_point_region[0],"\t",switch_point_region[$-1],"\t","-","\t",patient_switch_points,"\t",father_switch_points,"\t",mother_switch_points,"\t",perc_patient,"\t",perc_father,"\t",perc_mother,"\t","FAILED","\t","PATIENT PERCENTAGE","\t",patient_uniq,"\t",father_uniq,"\t",mother_uniq);
    } else if (patient_switch_points < patient_count) {
      output.writeln(chr,"\t",switch_point_region[0],"\t",switch_point_region[$-1],"\t","-","\t",patient_switch_points,"\t",father_switch_points,"\t",mother_switch_points,"\t",perc_patient,"\t",perc_father,"\t",perc_mother,"\t","FAILED","\t","PATIENT COUNT","\t",patient_uniq,"\t",father_uniq,"\t",mother_uniq);
    } else {
      double tuned = 0.0;
//       double tuned = fine_tuning(chr,switch_point_region[0], switch_point_region[$-1], samples[patientName].keys);
      output.writeln(chr,"\t",switch_point_region[0],"\t",switch_point_region[$-1],"\t",to!ulong(tuned),"\t",patient_switch_points,"\t",father_switch_points,"\t",mother_switch_points, "\t",perc_patient,"\t",perc_father,"\t",perc_mother,"\t","PASS","\t","-","\t",patient_uniq,"\t",father_uniq,"\t",mother_uniq);
    }
  }  
}

double fine_tuning( string chr, ulong start, ulong end, string[] files) {
//   writeln(chr,"\t",start,"\t",end);
  double[] tuned_start;
  double[] tuned_end;
  foreach( file ; files ) {
    File inFile = File(file,"r");
    while(!inFile.eof()) {
      string[] line = chomp(inFile.readln()).split("\t");
      if (line.length > 1 && line[0] != "#CHR") {
	if (line[0] == chr && to!ulong(line[1]) >= start && to!ulong(line[1]) <= end) {
	  ulong[2] tuned_positions = get_regions(file, chr, to!ulong(line[1]));
	  if (tuned_positions[0] != 0) {
	    tuned_start ~= to!double(tuned_positions[0]);
	  }
	  if (tuned_positions[1] != 0) {
	    tuned_end ~= to!double(tuned_positions[1]);
	  }
	}
      }
    }
  }
  double[] tuned = tuned_start ~ tuned_end;
  double tuned_mean = mean(tuned);
  return(tuned_mean);
//   writeln(chr,"\t",start,"\t",end,"\t",tuned_positions.keys.sort[0],"\t",tuned_positions.keys.sort[$-1]);
//   writeln(tuned_positions);
//   writeln();
}

ulong[2] get_regions( string file, string chr, ulong pos ) {
  file = file.replace("points","bins");
  File inFile = File(file,"r");
  string[] before_bin2;  
  string[] before_bin;
  string[] after_bin;
  string[] after_bin2;  
  while(!inFile.eof()) {
    string[] line = chomp(inFile.readln()).split("\t");
    if (line.length > 0 && line[0] == chr) {
      if (to!ulong(line[2]) <= pos) {
	if (before_bin2.length == 0) {
	  before_bin2 = [line[1],line[2],line[6]];
	} else if (before_bin.length == 0) {
	  before_bin = [line[1],line[2],line[6]];
	} else {
	  before_bin2 = before_bin;
	  before_bin = [line[1],line[2],line[6]];
	}
      }
      if (to!ulong(line[1]) >= pos) {
	if (after_bin.length == 0) {
	  after_bin = [line[1],line[2],line[6]];
	} else if (after_bin2.length == 0) {
	  after_bin2 = [line[1],line[2],line[6]];
	} else {
	  break;
	} 
      }
    }
  }
//   writeln(file,"\t",chr,"\t",pos,"\t",before_bin2,"\t",before_bin,"\t",after_bin,"\t",after_bin2);
//   writeln(before_bin[2],"\t",after_bin[2]);
  ulong[2] tuned_positions = get_tuned_position( file, chr, pos, to!uint(before_bin2[0]), to!uint(after_bin2[1]), before_bin[2], after_bin[2]);
  return(tuned_positions);
}

ulong[2] get_tuned_position( string file, string chr, ulong pos, uint start, uint end, string from, string to ) {
  ulong[2] tuned_positions;  
  file = file.replace(".bins","");
  BamReader bam = new BamReader(file);
  auto reads = bam[chr][start..end];
  BamRead[] subreads;
  foreach( read ; reads ) {
    if (read.mapping_quality > MQ) {
      if (subreads.length < 14) {
	subreads ~= read;
      } else {
	int[2] before_strand = [0, 0];
	int[2] after_strand = [0, 0];
	for (int i; i < subreads.length; i++) {
	  if (subreads[i].strand == '+' && i < 7) {
	    before_strand[0]++;
	  } else if (subreads[i].strand == '-' && i < 7) {
	    before_strand[1]++;
	  } else if (subreads[i].strand == '+' && i >= 7) {
	    after_strand[0]++;
	  } else if (subreads[i].strand == '-' && i >= 7) {
	    after_strand[1]++;
	  }
	}
	if (from == "F" && to == "R" && before_strand[0] == 7 && after_strand[1] == 7) {
	  tuned_positions = [subreads[6].position,subreads[7].position];
// 	  writeln(file,"\t",chr,"\t",pos,"\t",start,"\t",end,"\t",from,"\t",to,"\t",subreads[6].position,"\t",subreads[7].position);
	} else if (from == "F" && to == "FR" && before_strand[0] == 7 && after_strand[1] == 4) {
	  tuned_positions = [subreads[6].position,subreads[7].position];
// 	  writeln(file,"\t",chr,"\t",pos,"\t",start,"\t",end,"\t",from,"\t",to,"\t",subreads[6].position,"\t",subreads[7].position);
	} else if (from == "R" && to == "F" && before_strand[1] == 7 && after_strand[0] == 7) {
	  tuned_positions = [subreads[6].position,subreads[7].position];
// 	  writeln(file,"\t",chr,"\t",pos,"\t",start,"\t",end,"\t",from,"\t",to,"\t",subreads[6].position,"\t",subreads[7].position);
	} else if (from == "R" && to == "FR" && before_strand[1] == 7 && after_strand[0] == 4) {
	  tuned_positions = [subreads[6].position,subreads[7].position];
// 	  writeln(file,"\t",chr,"\t",pos,"\t",start,"\t",end,"\t",from,"\t",to,"\t",subreads[6].position,"\t",subreads[7].position);
	} else if (from == "FR" && to == "R" && before_strand[0] == 4 && after_strand[1] == 7) {
	  tuned_positions = [subreads[6].position,subreads[7].position];
// 	  writeln(file,"\t",chr,"\t",pos,"\t",start,"\t",end,"\t",from,"\t",to,"\t",subreads[6].position,"\t",subreads[7].position);
	} else if (from == "FR" && to == "F" && before_strand[1] == 4 && after_strand[0] == 7) {
	  tuned_positions = [subreads[6].position,subreads[7].position];
// 	  writeln(file,"\t",chr,"\t",pos,"\t",start,"\t",end,"\t",from,"\t",to,"\t",subreads[6].position,"\t",subreads[7].position);
	} 
	subreads = subreads[1..$-1];
	subreads ~= read;
      }
    }
  }
  return(tuned_positions);
}

string get_output_name( string folder ) {
  if (folder.endsWith("/")) {
    folder = replace(folder,regex(r"\/$"),"");
  }
  string[] path = folder.split("/");  
  string outName = path[$-1];
  return outName;
}

string process_slice_reads(BamReadBlock[] reads, string chr, ulong fwd_reads, ulong rev_reads) {
  int start = reads[0].position;
  int end = reads[$-1].position;
  double fwd_ratio = to!double(fwd_reads)/(to!double(rev_reads)+to!double(fwd_reads));
  string category = "FR";
  if (fwd_ratio > 0.99) {
    category = "F";
  } else if (fwd_ratio < 0.01) {
    category = "R";
  }
  return(category);
}

double mean(double[] array) {
  double sum = 0.0;
  foreach( item ; array ) {
    sum += item;
  }
  double mean = sum/to!double(array.length);
  return(mean);
}

double median(double[] array) {
  array.sort;
  double pos = to!double(array.length/2);
  double median = 0;
  if ((array.length % 2) != 0) {
    median = array[to!ulong(pos+0.5)];
  } else {
    median = (array[to!ulong(pos)]+array[to!ulong(pos-1)])/2;
  }
  return(median);
}

double stdev(double[] array) {
  double mean = mean(array);
  double sqtotal = 0.0;
  foreach( item ; array ) {
    sqtotal += (mean-to!double(item))^^2;
  }
  auto std = (sqtotal/array.length)^^0.5;
  return(std);
}

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

class Point {
  string chr;
  ulong pos;
  this(string chr, ulong pos) {
    this.chr = chr;
    this.pos = pos;
  }
}
