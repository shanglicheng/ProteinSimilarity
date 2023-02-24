#!/usr/bin/perl
my $matrix_index_temp = shift @ARGV;
my $matrix_index_tag = shift @ARGV;
my $pdb_name = shift @ARGV;
$_ = $pdb_name;
s/\.pdb//;
my $out_index = $matrix_index_temp."/".$matrix_index_tag."_character_index/";
if(-e $out_index){
} else {
	mkdir ($out_index);
}
for ( my $aa=0;$aa<5;$aa++ ) {
	my $matrix_index = $matrix_index_temp."/".$matrix_index_tag."_linmatrix_".$aa."/";
	undef @filenames;
	opendir(FILENAME,"$matrix_index");
	@filenames = readdir FILENAME;
	close(FILENAME);
	foreach  (@filenames) {
		if ( !($_ =~ /^\./) ) {
			$path = $matrix_index.$_;
			s/\.linjie\.matrix//;
			$out = $out_index.$_.".cindex";
			open(OUT,">>$out");
			open(FILE,"$path");
			@files = <FILE>;
			close(FILE);
			undef(@a_temp);
			undef($col_temp);
			for ( my $k=0;$k<@files;$k++ ) {
				$a_temp[$k] = $files[$k];
				chomp($a_temp[$k]);
				$a_temp[$k] = [(split/\s+/,$a_temp[$k])];
			}
			$col_temp = $files[0];
			undef(@a);
			chomp($col_temp);
			@col_temp_ = split/\s+/,$col_temp;
			$col = @col_temp_;
			$row = @files;
			undef($b_temp);
			for ( my $n=0;$n<$col;$n++ ) {
				for ( my $m=0;$m<$row;$m++ ) {
					$b_temp[$n][$m] = $a_temp[$m][$n];
				}
			}
			for ( my $n=0;$n<$row;$n++ ) {
				for ( my $ai=0;$ai<$row;$ai++) {
					$a[$n][$ai] = 0;
					for ( my $m=0;$m<$col;$m++ ) {
						$a[$n][$ai] = $a[$n][$ai] + $a_temp[$n][$m]*$b_temp[$m][$ai];
					}
				}
			}
			for ( my $kk=0;$kk<$row;$kk++ ) {
				$x0[$kk] = 1;
			}
			my $yi = 0;
			undef $x_temp;
			foreach  (@x0) {
				$x_temp = $_;
				$temp = &maxa(@x0);
				$y[$yi] = $x_temp/$temp;
				$yi++;
			}
			my $x1i = 0;
			for ( my $i=0;$i<@y ;$i++) {
				$x1[$x1i] = 0;
				for (my $j=0;$j<@y ;$j++ ) {
					$x1[$x1i] = $x1[$x1i] + $a[$i][$j]*$y[j];
				}
				$x1i++;
			}
			$max_x1 = &maxa(@x1);
			$max_x0 = &maxa(@x0);
			$max_x1_x0 = abs($max_x1-$max_x0);
			while ( $max_x1_x0>0.001 ) {
				for (my $i=0;$i<@x0 ;$i++ ) {
					$x0[$i] = $x1[$i];
				}
				my $yi = 0;
				foreach  (@x0) {
					$x_temp = $_;
					$temp = &maxa(@x0);
					$y[$yi] = $x_temp/$temp;
					$yi++;
				}
				my $x1i = 0;
				for ( my $i=0;$i<@y;$i++ ) {
					$x1[$x1i] = 0;
					for (my $j=0;$j<@y ;$j++ ) {
						$x1[$x1i] = $x1[$x1i] + $a[$i][$j]*$y[$j];
					}
					$x1i++;
				}
				$max_x1 = &maxa(@x1);
				$max_x0 = &maxa(@x0);
				$max_x1_x0 = abs($max_x1-$max_x0);
			}
			$eig = &maxa(@x1);
			syswrite OUT, $eig."\n";
		}
	}
}
sub maxa{
	@x = @_;
	$k = 0;
	for ( my $i=1;$i<@x;$i++ ) {
		$xi = abs($x[$i]);
		$xk = abs($x[$k]);
		if ( $xi > $xk ) {
			$k = $i;
		}
	}
	$y = $x[$k];
}