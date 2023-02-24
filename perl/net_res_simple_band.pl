#!/usr/bin/perl
my $vanresnet_index = shift @ARGV;
my $out_index = $vanresnet_index;
my $pdb_name = shift @ARGV;
$_ = $pdb_name;
s/\.pdb//;
@files = $_.".vanresnet";
foreach  (@files) {
	if ( !($_ =~ /^\./) ) {
		$filepath = $vanresnet_index."/".$_;
		s/\.vanresnet//;
		$outpath = $out_index."/".$_.".simplenet";
		open(OUT,">$outpath");
		open(FILE,"$filepath");
		@atoms = <FILE>;
		close(FILE);
		print @atoms."\n";
		for ( my $i=0;$i<@atoms;$i++ ) {
			@i_temp_ = split/\s+/,$atoms[$i];
			@res_0_i = split/\_/,$i_temp_[0];
			@res_2_i = split/\_/,$i_temp_[2];
			@atom_i_ = split/\_/,$i_temp_[1];
			$aa_tag = 0;#
			for ( my $j=$i;$j<@atoms+1;$j++ ) {
				@j_temp_ = split/\s+/,$atoms[$j];
				@res_0_j = split/\_/,$j_temp_[0];
				@res_2_j = split/\_/,$j_temp_[2];
				@atom_j_ = split/\_/,$j_temp_[1];
				if ( $i_temp_[0] ne $j_temp_[0] ) {
					$i = $j-1;
					last;
				} else {
					if ( ($aa_tag < 1)&&(abs($res_0_j[2]-$res_2_j[2])==1) ) {
						syswrite OUT, $j_temp_[0]." ".$j_temp_[3]." ".$j_temp_[2]."\n";#" ".$atom_j_[0]."_".$atom_j_[3].
						$aa_tag++;
					}
				}
			}
		}
		syswrite OUT, "\n";
		for ( my $i=0;$i<@atoms;$i++ ) {
			@i_temp_ = split/\s+/,$atoms[$i];
			@res_0_i = split/\_/,$i_temp_[0];
			@res_2_i = split/\_/,$i_temp_[2];
			@atom_i_ = split/\_/,$i_temp_[1];
			$rr_num = 0;#
			$i_temp = $i;
			undef(@rr_res_temp);
			for ( my $j=$i;$j<@atoms+1;$j++ ) {
				@j_temp_ = split/\s+/,$atoms[$j];
				@res_0_j = split/\_/,$j_temp_[0];
				@res_2_j = split/\_/,$j_temp_[2];
				@atom_j_ = split/\_/,$j_temp_[1];
				if ( $i_temp_[0] ne $j_temp_[0] ) {
					$i = $j-1;
					last;
				} else {
					if ( abs($res_0_j[2] - $res_2_j[2]) > 1 ) {
						$rr_res[$rr_num] = $j_temp_[2];
						$rr_num++;
					}
				}
			}
			for ( my $r=0;$r<$rr_num;$r++) {
				$rr_res_temp[$r] = $rr_res[$r];
			}
			$i_end = $i;
			my %saw;
			@saw{@rr_res_temp} = ();
			my @uniq_res = sort keys %saw;
			for ( my $j=$i_temp;$j<$i_end;$j++ ) {
				@j_temp_ = split/\s+/,$atoms[$j];
				@res_0_j = split/\_/,$j_temp_[0];
				@res_2_j = split/\_/,$j_temp_[2];
				@atom_j_ = split/\_/,$j_temp_[1];
				for ( my $rr=0;$rr<@uniq_res;$rr++ ) {
					if ( $j_temp_[2] eq $uniq_res[$rr] ) {
						syswrite OUT, $j_temp_[0]." ".$j_temp_[3]." ".$j_temp_[2]."\n";#" ".$atom_j_[0]."_".$atom_j_[3].
						$uniq_res[$rr] = "";
					}
				}
			}
		}
	}
}