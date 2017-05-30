#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$FindBin::Bin/../source";
use CF::Constants;
use CF::Helpers;
use File::Copy qw(move);

##########################################################################
# Copyright 2014, Philip Ewels (phil.ewels@scilifelab.se)                #
# Copyright 2016, Russell S, Hamilton (rsh46@cam.ac.uk) (doublemap)      #
#                                                                        #
# This file is part of Cluster Flow.                                     #
#                                                                        #
# Cluster Flow is free software: you can redistribute it and/or modify   #
# it under the terms of the GNU General Public License as published by   #
# the Free Software Foundation, either version 3 of the License, or      #
# (at your option) any later version.                                    #
#                                                                        #
# Cluster Flow is distributed in the hope that it will be useful,        #
# but WITHOUT ANY WARRANTY; without even the implied warranty of         #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          #
# GNU General Public License for more details.                           #
#                                                                        #
# You should have received a copy of the GNU General Public License      #
# along with Cluster Flow.  If not, see <http://www.gnu.org/licenses/>.  #
##########################################################################

# Module requirements
my %requirements = (
	'cores' 	=> ['1', '10'],
	'memory' 	=> ['8G', '15G'],
	'modules' 	=> ['tophat', 'samtools'],
	'references'=> 'bowtie2',
	'time' 		=> sub {
		my $cf = $_[0];
		my $num_files = $cf->{'num_starting_merged_aligned_files'};
		$num_files = ($num_files > 0) ? $num_files : 1;
		# Alignment can take ages. Be conservative..
		#return CF::Helpers::minutes_to_timestamp ($num_files * 10 * 60);
                return CF::Helpers::minutes_to_timestamp ($num_files * 90 * 60);
	}
);

# Help text
my $helptext = "".("-"x15)."\n Tophat Module\n".("-"x15)."\n
TopHat is a fast splice junction mapper for RNA-Seq reads.
It can run with or without a supplied GTF path but requires
a bowtie2 path.

This module fixes a bug in tophat that relates to MAPQ values;
by using –g 1 (as in the previous tophat module), ALL alignments will
get a MAPQ score of 50, irrespective of whether they were unique
alignments or not. This module uses –g 2 instead, so reporting up
to 2 alignments for a given read. One of the alignments is considered
the primary alignment, and the other one is flagged as secondary
alignment, and both alignments now get a MAPQ value of 3 to indicate
that there were multiple alignments. After Tophat has finished, this
module uses samtools to filter out the secondary alignments.

In short, this module produces a BAM file with an equivalent number
of hits as the old module, but with fixed MAPQ values.\n\n";

# Setup
my %cf = CF::Helpers::module_start(\%requirements, $helptext);

# MODULE
# Check that we have a genome defined
if(!defined($cf{'refs'}{'bowtie2'})){
	die "\n\n###CF Error: No bowtie2 path found in run file $cf{run_fn} for job $cf{job_id}. Exiting..";
} else {
	warn "\nAligning against bowtie2 path: $cf{refs}{bowtie2}\n\n";
}

# Use a GTF file if we have one
my $gtf = '';
if(defined($cf{'refs'}{'gtf'}) && -e $cf{'refs'}{'gtf'}){
	$gtf = "-G $cf{refs}{gtf}";
} elsif(defined($cf{'refs'}{'gtf'})){
	warn "\nWarning! GTF reference ".$cf{'refs'}{'gtf'}." specified, but file doesn't exist. Ignoring.\n\n";
}

open (RUN,'>>',$cf{'run_fn'}) or die "###CF Error: Can't write to $cf{run_fn}: $!";

# Print version information about the module.
warn "---------- Tophat version information ----------\n";
warn `tophat --version`;
warn "\n------- End of Tophat version information ------\n";

# Separate file names into single end and paired end
my ($se_files, $pe_files) = CF::Helpers::is_paired_end(\%cf, @{$cf{'prev_job_files'}});

# FastQ encoding type. Once found on one file will assume all others are the same
my $encoding = 0;

# Go through each single end files and run Tophat
my $output_dir;
if($se_files && scalar(@$se_files) > 0){
	foreach my $file (@$se_files){

		my $timestart = time;

		# Figure out the encoding if we don't already know
		if(!$encoding){
			($encoding) = CF::Helpers::fastq_encoding_type($file);
		}
		# Tophat only accepts --solexa-quals, --solexa1.3-quals and --integer-quals
		# Defaults to phred-33 so leave blank if that returned
		my $enc = "";
		if($encoding eq 'solexa'){
			$enc = '--solexa-quals';
		} elsif($encoding eq 'phred64'){
			$enc = '--solexa1.3-quals';
		} elsif($encoding eq 'integer'){
			$enc = '--integer-quals';
		}

		$output_dir = $file;
		$output_dir =~ s/.gz$//;
		$output_dir =~ s/.fq$//;
		$output_dir =~ s/.fastq$//;
		$output_dir =~ s/_[1-4]$//;
		$output_dir =~ s/_R[1-4]//;
		$output_dir .= '_tophat';

                my $output_dir2 = $output_dir . "/tophat_secondrun"; 
		my $output_fn   = $output_dir2."/accepted_hits_merged.bam";

		# Changing the command to -g 2 so that we can get proper MAPQ values
		$gtf = " -G $cf{refs}{bowtie2}.gtf";

		my $command1 = "tophat -p $cf{cores} -g 2 $enc $gtf -o $output_dir $cf{refs}{bowtie2} $file";
                #warn "\n###CFCMD $command1\n\n";
                $command1 = "echo '\n###CFCMD $command1\n\n 1>&2'; $command1";

                my $command2 = "bed_to_juncs < $output_dir/junctions.bed | sort -k 1,4 -u | sort -k 1,1 > $output_dir/junctions.junc";
                #warn "\n###CFCMD $command2\n\n";
                $command2 = "echo '\n###CFCMD $command2\n\n 1>&2'; $command2";

                my $command3 = "samtools fastq $output_dir/unmapped.bam | gzip -c > $output_dir/unmapped.fq.gz";
                #warn "\n###CFCMD $command3\n\n";
                $command3 = "echo '\n###CFCMD $command3\n\n 1>&2'; $command3";

                my $command4 = "tophat -p $cf{cores} -g 2 $enc $gtf -o $output_dir2 -j $output_dir/junctions.junc --no-novel-junc $cf{refs}{bowtie2} $output_dir/unmapped.fq.gz";
                #warn "\n###CFCMD $command4\n\n";
                $command4 = "echo '\n###CFCMD $command4\n\n 1>&2'; $command4";

                my $command5 = "samtools merge $output_fn $output_dir/accepted_hits.bam $output_dir2/accepted_hits.bam";
		#warn "\n###CFCMD $command5\n\n";
                $command5 = "echo '\n###CFCMD $command5\n\n 1>&2'; $command5";

		if( !system($command1) and !system($command2) and !system($command3) and !system($command4) and !system($command5)){
			# Tophat worked - print out resulting filenames
			my $duration =  CF::Helpers::parse_seconds(time - $timestart);
			warn "###CF Tophat (SE mode) successfully exited, took $duration..\n";
			if(-e $output_fn){
                                $output_fn = clean_output($output_dir, $output_dir2, "single");
				print RUN "$cf{job_id}\t$output_fn\n";
			} else {
				warn "\n###CF Error! Tophat output file $output_fn not found..\n";
			}
		} else {
			warn "\n###CF Error! Tophat (SE mode) failed for input file '$file': $? $!\n\n";
		}

	}
}

# Go through the paired end files and run Tophat
if($pe_files && scalar(@$pe_files) > 0){
	foreach my $files_ref (@$pe_files){
		my @files = @$files_ref;
		if(scalar(@files) == 2){

			my $timestart = time;

			# Figure out the encoding if we don't already know
			if(!$encoding){
				($encoding) = CF::Helpers::fastq_encoding_type($files[0]);
			}
			# Tophat only accepts --solexa-quals, --solexa1.3-quals and --integer-quals
			# Defaults to phred-33 so leave blank if that returned
			my $enc = "";
			if($encoding eq 'solexa'){
				$enc = '--solexa-quals';
			} elsif($encoding eq 'phred64'){
				$enc = '--solexa1.3-quals';
			} elsif($encoding eq 'integer'){
				$enc = '--integer-quals';
			}

			$output_dir = $files[0];
			$output_dir =~ s/.gz$//;
			$output_dir =~ s/.fq$//;
			$output_dir =~ s/.fastq$//;
			$output_dir =~ s/_[1-4]$//;
			$output_dir =~ s/_R[1-4]//;
			$output_dir .= '_tophat';

                        my $output_dir2 = $output_dir . "/tophat_secondrun";

			#my $output_fn = $output_dir."/accepted_hits.bam";
                        my $output_fn = $output_dir2."/accepted_hits_merged.bam";

                        # Changing the command to -g 2 so that we can get proper MAPQ values
                        $gtf = " -G $cf{refs}{bowtie2}.gtf";

			my $command1 = "tophat -p $cf{cores} -g 2 $enc $gtf -o $output_dir $cf{refs}{bowtie2} $files[0] $files[1]";
			warn "\n###CFCMD $command1\n\n";
                        $command1 = "echo '\n###CFCMD $command1\n\n'; $command1";

                        my $command2 = "bed_to_juncs < $output_dir/junctions.bed | sort -k 1,4 -u | sort -k 1,1 > $output_dir/junctions.junc";
                        #warn "\n###CFCMD $command2\n\n";
                        $command2 = "echo '\n###CFCMD $command2\n\n'; $command2";

                        my $command3 = "samtools fastq -1 $output_dir/unmapped_1.fq -2 $output_dir/unmapped_2.fq $output_dir/unmapped.bam; " .
                                       "gzip $output_dir/unmapped_1.fq $output_dir/unmapped_2.fq";
                        #warn "\n###CFCMD $command3\n\n";
                        $command3 = "echo '\n###CFCMD $command3\n\n'; $command3";

                        my $command4 = "tophat -p $cf{cores} -g 2 $enc $gtf -o $output_dir2 -j $output_dir/junctions.junc --no-novel-junc $cf{refs}{bowtie2} " . 
                                       "$output_dir/unmapped_1.fq.gz $output_dir/unmapped_2.fq.gz";
                        #warn "\n###CFCMD $command4\n\n";
                        $command4 = "echo '\n###CFCMD $command4\n\n'; $command4";

                        my $command5 = "samtools merge $output_fn $output_dir/accepted_hits.bam $output_dir2/accepted_hits.bam";
                        #warn "\n###CFCMD $command5\n\n";
                        $command5 = "echo '\n###CFCMD $command5\n\n'; $command5";

                        if( !system($command1) and !system($command2) and !system($command3) and !system($command4) and !system($command5)){
				# Tophat worked - print out resulting filenames
				my $duration =  CF::Helpers::parse_seconds(time - $timestart);
				warn "###CF Tophat (PE mode) successfully exited, took $duration....\n";
				if(-e $output_fn){
					$output_fn = clean_output($output_dir, $output_dir2, "paired");
					print RUN "$cf{job_id}\t$output_fn\n";
				} else {
					warn "\n###CF Error! Tophat output file $output_fn not found\n";
				}
			} else {
				warn "\n###CF Error! Tophat (PE mode) failed for input file '$files[0]': $? $!\n\n";
			}
		} else {
			warn "\n###CF Error! Tophat paired end files had ".scalar(@files)." input files instead of 2\n";
		}

	}
}


# Clear up Tophat output
sub clean_output {

    my ($output_dir, $output_dir2, $readstate) = @_;
    # e.g.
    # 1. SLX-9615.A010.C66Y8ANXX.s_5.r_1_trimmed_tophat
    # 2. /tophat_secondrun
    # Main results in dir2
    # accepted_hits_merged.bam

    my $results_fn     = $output_dir.".accepted_hits_merged.bam";
    my $results_fn_sum = $output_dir.".align_summary.txt";

    if($readstate eq "paired")
      {
        $results_fn     = $output_dir.".accepted_hits_merged_TMP.bam";
      }

    if(-e "$output_dir2/accepted_hits_merged.bam" && -e "$output_dir/align_summary.txt"){

        # Use samtools to remove secondary alignments
    	my $samtools_command  = "samtools view -F 0x100 -b -h $output_dir2/accepted_hits_merged.bam | samtools sort -n > $results_fn";

        warn "Using samtools to filter out non-primary alignments from the tophat alignments\n\n";
        warn "\n###CFCMD $samtools_command\n\n";
    	if(system ($samtools_command))
          {
    	    $results_fn = "";
            warn "###CF Error! Failed to remove non-primary alignments from $output_dir2/accepted_hits_merged.bam to ${output_dir2}.accepted_hits_merged.bam: $? $!\n";
      	  } 
        else 
         {
            # System returned 0, samtools worked..
            warn "Successfully created new alignments file without duplicates.\n"; 
         }


        my $results_fn_fixed;
        my $fixmate_command;
        my $viewf0x01_command;
        if($readstate eq "paired")
          {
            $results_fn_fixed  = $output_dir.".accepted_hits_merged_TMP_fixed.bam";
            $fixmate_command   = "samtools fixmate -r $results_fn $results_fn_fixed";

            warn "\n###CFCMD $fixmate_command\n\n";

            if(system($fixmate_command))
              {
                $results_fn = "";
                warn "###CF Error! PE fixmate failed for $results_fn > $results_fn_fixed: $? $!\n"; 
              }
            else 
              {
                unlink($results_fn); 
                warn "Successfully ran fixmate on $results_fn\n";
              }

            $results_fn        = $output_dir.".accepted_hits_merged.bam";
            $viewf0x01_command = "samtools view -f 0x01 -b $results_fn_fixed > $results_fn";
            warn "\n###CFCMD $viewf0x01_command\n\n";
 
            if(system($viewf0x01_command))
              {
                $results_fn = "";
                warn "###CF Error! PE viewf0x01 failed for $results_fn_fixed: $? $!\n";
              }
            else 
             { 
                unlink($results_fn_fixed);
                warn "Successfully ran viewf0x01 on $results_fn\n";
             }


          }


            # Delete old alignments file with the duplicates
    #        if(unlink "$output_dir/accepted_hits.bam"){
    #            warn "Deleted original Tophat alignment file with duplicates.\n"
    #        } else {
    #            warn "Could not delete original Tophat alignment file with duplicates.\n"
    #        }

# First run summary file
if(move ("$output_dir/align_summary.txt", $output_dir.".tophatrun1.align_summary.txt")){
    warn "Moved tophat summary to ${output_dir}.tophatrun1.align_summary.txt: $!\n";
} else  {
    warn "Could not move tophat run 1 summary to ${output_dir}.tophatrun1.align_summary.txt: $!\n"; }
# Second run summary file
if(move ("$output_dir2/align_summary.txt", $output_dir.".tophatrun2.align_summary.txt")){
    warn "Moved tophat summary to ${output_dir}.tophatrun2.align_summary.txt: $!\n";
} else  {
    warn "Could not move tophat run 2 summary to ${output_dir}.tophatrun2.align_summary.txt: $!\n"; }

# REMOVE THE UNMEREGED BAM FILES

if(unlink glob "$output_dir/logs/*"){
  warn "Deleted other tophat output files: $output_dir/logs/*\n";
} else {
  warn "Could not delete other tophat log files: $output_dir/logs/*: $!\n";}
if(rmdir "$output_dir/logs"){
  warn "Deleted tophat logs directory $output_dir/logs\n";
} else {
  warn "Could not delete tophat logs directory $output_dir/logs: $!\n"; }

if(unlink glob "$output_dir2/logs/*"){
  warn "Deleted other tophat output files: $output_dir2/logs/*\n";
} else {
  warn "Could not delete other tophat log files: $output_dir2/logs/*: $!\n";}
if(rmdir "$output_dir2/logs"){
  warn "Deleted tophat logs directory $output_dir2/logs\n";
} else {
  warn "Could not delete tophat logs directory $output_dir2/logs: $!\n"; }


            # Go through and clear up all of the other tophat files
            #if(move ("$output_dir/align_summary.txt", $output_dir.".align_summary.txt")){
            #    warn "Moved filtered tophat summary to ${output_dir}.align_summary.txt\n";
            #    if(unlink glob "$output_dir/logs/*"){
            #        warn "Deleted other tophat output files: $output_dir/logs/*\n";
            #        if(rmdir "$output_dir/logs"){
            #            warn "Deleted tophat logs directory $output_dir/logs\n";
 #                       if(unlink glob "$output_dir/*"){
 #                           warn "Deleted other tophat output files: $output_dir/*\n";
 #                           if(rmdir $output_dir){
 #                               warn "Deleted tophat output directory $output_dir\n";
 #                           } else {
 #                               warn "Could not delete tophat output directory $output_dir: $!\n";
 #                           }
 #                       } else {
 #                           warn "Could not delete other tophat output files: $output_dir/*: $!\n";
 #                       }
            #        } else {
            #            warn "Could not delete tophat logs directory $output_dir/logs: $!\n";
            #        }
            #    } else {
            #        warn "Could not delete other tophat log files: $output_dir/logs/*: $!\n";
            #    }
            #} else {
            #    warn "Could not move tophat summary to ${output_dir}.align_summary.txt: $!\n";
            #}
  #      }
  #  } else {
  #      warn "Could not find tophat output files $output_dir/accepted_hits.bam and $output_dir/align_summary.txt\n";
    }

    return $results_fn;
}


close (RUN);
