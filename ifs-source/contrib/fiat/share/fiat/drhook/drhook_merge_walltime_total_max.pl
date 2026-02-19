#
# drhook_merge_walltime.pl
#
# For merging wall clock time results from different MPI-tasks
# i.e. DR_HOOK_OPT=prof 
#
# Original script by Eckhard Tschirschnitz, Cray, 2006 (Mflop/s)
#
# Usage: cat drhook.* | perl -w drhook_merge_walltime_total_max.pl
#
# (sorts w.r.t max-value, not average)
#

use strict;

# this expects concatenated dr_hook listings (wall clock time listings)

# note below: self in fact denotes self + children times i.e. total time

my $bignum = 999999999;

my $skip = 1;
#my $threshold = 1.0;
my $threshold = 0.001;
my $tottim = 0;

my $maxwall = 0;
my $minwall = $bignum;
my $avgwall = 0;
my $stdevwall = 0;

my $nproc = 0; # no of MPI-tasks
my $omp = 0; # max. no of OpenMP-threads encountered
my $exe = ""; # the name of the executable

my %namelist = ();

my %sumself = ();
my %sum2self = ();
my %maxself = ();
my %minself = ();
my %ompself = ();

my %numcalls = ();

for (<>) {
  chomp; # get rid of the newline character

  next if (m/^\s*$/); # a blank line

  if (m/^\s*Profiling\b/) { 
    if ($nproc == 0) {
      $exe = $1 if (m/program='([^\'].*)'/);
    }
    $nproc++;
    $skip = 1; 
    next; # for (<>)
  }
  elsif (m/^\s+Wall-time\s+is\s+(\S+)\s+/) {
    my $value = $1;
    $maxwall = $value if ($value > $maxwall);
    $minwall = $value if ($value < $minwall);
    $avgwall += $value;
    $stdevwall += $value * $value;
    next; # for (<>)
  }
  elsif (m/^\s+1\s+/) {
    foreach my $name (keys %ompself) {
      my $self = $ompself{$name};
      $sumself{$name} += $self;
      $sum2self{$name} += $self * $self;
      $maxself{$name} = $self if ($self > $maxself{$name});
      $minself{$name} = $self if ($self < $minself{$name});
      $tottim += $self;
    }
    %ompself = ();
    $skip = 0;
  }

  if ($skip == 0) {
    #         rank  %time cumul self   total   #_of_calls self:ms/call tot:ms/call routine_name
    if (m/^\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+\S+\s+\S+\s+(.*)/) {
      my $self = $1;
      my $ncalls = $2;
      my $name = $3;
      my $tid = 0;
      $name =~ s/\s+//g;
      $name =~ s/^[*]//;
      #print "$self $name\n";
      if ($name =~ m/^(.*)[@](\d+)/) {
	$tid = $2;
	$omp = $tid if ($tid > $omp);
	$name = $1;
      }
      $namelist{$name} = $name if (!defined($namelist{$name}));
      if (!defined($sumself{$name})) {
	$sumself{$name} = 0;
	$sum2self{$name} = 0;
	$maxself{$name} = 0;
	$minself{$name} = $bignum;
	$numcalls{$name} = 0.0;
      }
      $numcalls{$name} += $ncalls;
      # Account the most expensive OpenMP thread only
      $ompself{$name} = 0 if (!defined($ompself{$name}));
      $ompself{$name} = $self if ($self > $ompself{$name});
    }
  }
}

if ($nproc > 0) {
  # One final time ...
  foreach my $name (keys %ompself) {
    my $self = $ompself{$name};
    $sumself{$name} += $self;
    $sum2self{$name} += $self * $self;
    $maxself{$name} = $self if ($self > $maxself{$name});
    $minself{$name} = $self if ($self < $minself{$name});
    $tottim += $self;
  }

  print STDOUT "Name of the executable : $exe\n";
  print STDOUT "Number of MPI-tasks : $nproc\n";
  print STDOUT  "Number of OpenMP-threads : $omp\n";
#  printf ("Total time : %.3f secs (sum over all MPI-tasks; only max OpenMP-time accounted for)\n", $tottim);
  $avgwall /= $nproc;
  if ($nproc > 1) {
    $stdevwall = ($stdevwall - $nproc * $avgwall * $avgwall)/($nproc - 1);
    $stdevwall = ($stdevwall > 0) ? sqrt($stdevwall) : 0; # be careful with rounding of errors
  }
  else {
    $stdevwall = 0;
  }
  printf STDOUT ("Wall-times over all MPI-tasks (secs) : Min=%.3f, Max=%.3f, Avg=%.3f, StDev=%.3f\n", 
		 $minwall, $maxwall, $avgwall, $stdevwall);
  printf ("Routines whose total time (i.e. sum) > %.3f secs will be included in the listing\n",$threshold);
  printf STDOUT ("%7s %10s %10s %10s %8s %8s %12s : %s\n",
		 "Avg-%", "Avg.time", 
		 "Min.time", "Max.time", "St.dev", 
		 "Imbal-%", 
		 "# of calls", "Name of the routine");
  # open(PIPE,"|sort -nr|sed 's/ %/%/'"); # sorts w.r.t average percentage (commented ayt)
  open(PIPE,"|sort -nr +4|sed 's/ %/%/'"); # sorts w.r.t. max column# (column number starts from 0, not 1 in unix sort)

  # Values accounted for
  my $acc_avgpercent = 0;
  my $acc_avgtime = 0;
  my $acc_maxtime = 0;
  my $acc_mintime = $bignum;

  foreach my $name (keys %sumself) {
    my $value = $sumself{$name};
    if ($value > $threshold) {
      my $percent = $value/$tottim*100;
      my $aveself = $value/$nproc;
      my $stdev = 0;
      if ($nproc > 1) {
	$stdev = $sum2self{$name};
	$stdev = ($stdev - $nproc * $aveself * $aveself)/($nproc - 1);
	$stdev = ($stdev > 0) ? sqrt($stdev) : 0; # be careful with rounding of errors
      }
      printf PIPE ("%6.2f %% %10.3f %10.3f %10.3f %8.3f %7.2f%% %12.0f : %s\n",
		   $percent, $aveself,
		   $minself{$name},$maxself{$name},$stdev,
		   ($maxself{$name} - $minself{$name})/$maxself{$name}*100,
		   $numcalls{$name}, $name);
      # Update values accounted for
      $acc_avgpercent += $percent;
      $acc_avgtime += $aveself;
      $acc_mintime = $minself{$name} if ($acc_mintime > $minself{$name});
      $acc_maxtime = $maxself{$name} if ($acc_maxtime < $maxself{$name});
    }
  }
  close(PIPE);
  printf STDOUT ("%6.2f%% %10.3f %10.3f %10.3f\n",
		 $acc_avgpercent,$acc_avgtime,
		 $acc_mintime, $acc_maxtime);
}

# A sample of typical output
__DATA__
Name of the executable : /fdb/mpm/RAPS9/fdb/eqgp/bin/ifsMASTER
Number of MPI-tasks : 64
Number of OpenMP-threads : 4
Wall-times over all MPI-tasks (secs) : Min=990.460, Max=1009.730, Avg=997.713, StDev=3.809
Routines whose total time (i.e. sum) > 1.000 secs will be included in the listing
  Avg-%   Avg.time   Min.time   Max.time   St.dev  Imbal-%   # of calls : Name of the routine
  4.63%     46.190      2.728     47.549    5.537   94.26%         1664 : >MPL-SCATTER_CTLVEC(524)
  3.22%     32.117     31.618     32.832    0.300    3.70%      1253376 : LWPTL
  3.10%     30.895     30.554     31.394    0.199    2.68%      1253376 : SWNIAD
  3.04%     30.287     29.697     31.075    0.287    4.43%      1253376 : LWPAD
  2.04%     20.350      7.867     21.129    1.933   62.77%      3034240 : BROADCREAL
  1.73%     17.288     17.098     17.562    0.088    2.64%    753223680 : SWDE
  1.73%     17.262     10.546     27.521    2.945   61.68%       158400 : >MPL-SLCOMM1_COMMS(509)
  1.61%     16.103     11.131     22.163    3.132   49.78%       316544 : >MPL-SLCOMM2A_COMMS(512)
  1.52%     15.172     14.947     15.551    0.120    3.88%      1253376 : SWNITL
  1.45%     14.488      0.024     16.144    2.017   99.85%       140544 : >MPL-IRCVGPF(527)
  1.45%     14.428     13.744     15.056    0.307    8.71%      2506752 : LWAI
  1.44%     14.383     13.958     15.184    0.297    8.07%    300810240 : SWDEAD
  1.43%     14.317     14.163     14.502    0.068    2.34%    300810240 : SWDETL
  1.36%     13.553     13.129     14.112    0.257    6.97%      6266880 : LAITRITLAD
  1.30%     12.980      6.884     15.539    1.590   55.70%    922794098 : CUADJTQ
  1.29%     12.909     12.566     13.452    0.178    6.59%      2506752 : SWRAD
