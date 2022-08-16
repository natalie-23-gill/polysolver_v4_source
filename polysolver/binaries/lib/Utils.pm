
=head1 LICENSE

Strelka Workflow Software
Copyright (c) 2009-2013 Illumina, Inc.

This software is provided under the terms and conditions of the
Illumina Open Source Software License 1.

You should have received a copy of the Illumina Open Source
Software License 1 along with this program. If not, see
<https://github.com/downloads/sequencing/licenses/>.

=cut


package Utils;

use base 'Exporter';

our @EXPORT = qw(
        errorX logX executeCmd checkFile checkDir checkMove
        getAbsPath checkMakeDir getBinList
        parseConfigIni writeConfigIni
    );

use warnings FATAL => 'all';
use strict;

use Carp;
use Cwd qw(realpath);
use File::Copy qw(move);
use File::Path qw(mkpath);


sub errorX($) {
    confess "\nERROR: " . $_[0] . "\n\n";
}

sub logX($) {
    print STDERR "INFO: " . $_[0] . "\n";
}


sub executeCmd($;$) {
    my $cmd = shift;
    my $isVerbose = shift;

    logX("Running: '$cmd'") if(defined($isVerbose) and $isVerbose);
    system($cmd) == 0
      or errorX("Failed system call: '$cmd'");
}

sub checkFile($;$) {
    my $file = shift;
    return if(-f $file);
    my $label = shift;
    errorX("Can't find" . (defined($label) ? " $label" : "") . " file: '$file'");
}

sub checkDir($;$) {
    my $dir = shift;
    return if(-d $dir);
    my $label = shift;
    errorX("Can't find" . (defined($label) ? " $label" : "") . " directory: '$dir'");
}

sub checkMove($$) {
    my ($old,$new) = @_;
    move($old,$new) || errorX("File move failed: $!\n\tAttempting to move '$old' to '$new'");
}




=item getAbsPath($path)

This procedure attempts to convert a path provided by the user on the
command line into an absolute path. It should be able to handle "~"
paths and conventional relative paths using ".." or ".". Resolution of
links should follow the convention of "Cwd::realpath".

B<Parameters:>

    $dirRef         - path (converted to absolute path in place)

B<Returns:>

    returns zero if successful, non-zero otherwise.

=cut
sub getAbsPath(\$) {
    my ($dirRef) = @_;
    my @tmp=glob($$dirRef);
    return 1 if(scalar(@tmp) != 1);
    my $ret = Cwd::realpath($tmp[0]);
    return 1 if !$ret && !($ret = File::Spec->rel2abs($tmp[0]));
    $$dirRef = $ret;
    return 0;
}



sub checkMakeDir($) {
    my $dir = shift;
    unless (-e $dir) {
        File::Path::mkpath($dir) || errorX("Can't create directory '$dir'");
    } else {
        errorX("Path is not a directory '$dir'\n") unless(-d $dir);
    }
}



sub getBinList($$) {
    my ($chromSize,$binSize) = @_;

    my $nm1 = (($chromSize-1) / $binSize);
    return [ map {sprintf("%04i",$_)} (0..$nm1) ];
}



sub parseConfigError($$) {
    my ($file,$line) = @_;
    errorX("Config file '$file' contains unexpected line '$line'\n");
}


sub parseConfigIni($) {
    my $file = shift;
    my %config;
    open(my $FH,"< $file") || errorX("Can't open config file '$file'");
    my $section = "noSection";
    while(<$FH>) {
        next if(/^[;#]/);
        next if(/^\s*$/);
        chomp;
        my $line=$_;
        my @ncl = split(/[;#]/);
        next unless(scalar(@ncl));
        my $nc = $ncl[0];
        if($nc =~ /^\s*\[([^\]]*)\]\s*$/) {
            $section = $1;
            next;
        }
        my ($key,$val) = map { s/^\s+//; s/\s+$//; $_ } split(/=/,$nc,2);
        unless(defined($key) && defined($val) && ($key ne "")) { parseConfigError($file,$line); }

        $config{$section}{$key} = $val;
    }
    close($FH);
    return \%config;
}


# minimal ini stringifier:
#
sub writeConfigIni($) {
    my $config = shift;
    my $val = "";
    for my $section (sort(keys(%$config))) {
        $val .= "\n[$section]\n";
        for my $key (sort(keys(%{$config->{$section}}))) {
            $val .= "$key = " . $config->{$section}{$key} . "\n";
        }
    }
    return $val;
}


1;
