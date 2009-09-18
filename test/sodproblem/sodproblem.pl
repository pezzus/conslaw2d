#!/usr/bin/perl

#main program
#controllo se gli argomenti sono corretti
$program = "sodproblem";
$numarg = @ARGV;
if (@ARGV < 1) {
	print("Usage: $program.pl [options] meshfile.msh\n");
	print("Options:\n");
	print("  --gnuplot\t\tGenerate plot and animation from gnuplot\n");
	print("  --interpolated\t\tInterpolate solution on vertices\n");
	exit();
}
my $meshfile;
$gnuplot = false;
$interpolated = false;
foreach $arg (@ARGV) {
	if ($arg eq "--gnuplot") {
		$gnuplot = true;
	} elsif($arg eq "--interpolated") {
		$interpolated = true;
	} else {
		$meshfile = $arg;
	}
}
#rimuovo dati di simulazioni precedenti
print ">> Rimozione dati precedenti ... \n";
system("rm data/* &> /dev/null");
system("rm output/* &> /dev/null");
#avvio la simulazione
print ">> Simulazione ... \n";
if(system("./$program ".join(" ",@ARGV))) {
	exit();
}
#generazione grafici
if ($gnuplot eq false) {
	exit();
}
@files = <data/solution*.dat>;
my $tmpname, $count = 0;
foreach $file (@files) {
$tmpname = substr($file,13,-4);
print ">> Genero grafico $count ...\n";
open(GNUPLOT,"| gnuplot");
print GNUPLOT <<gnuplot_commands;
	set term png transparent enhanced font "Vera,12" crop
	set output "output/$program-velocity-$tmpname.png"
	unset key
	unset xtics
	unset ytics
	set pm3d map
	set cbrange[0:2]
	splot '$file' using 1:2:(sqrt(\$4*\$4+\$5*\$5))
	set output "output/$program-density-$tmpname.png"
	set cbrange[0:2]
	splot '$file' using 1:2:3
gnuplot_commands
close GNUPLOT;
$count = $count + 1;
}
#animazione
print ">> Genero animazioni ... \n";
system ("convert -delay ".(0.05)." -loop 0 output/*density-*.png output/$program-density-movie.gif");
system ("convert -delay ".(0.05)." -loop 0 output/*velocity-*.png output/$program-velocity-movie.gif");
