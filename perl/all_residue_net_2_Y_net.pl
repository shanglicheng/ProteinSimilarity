#!/usr/bin/perl
my $simple_net_index = shift @ARGV;
my $outindex_temp_tag = shift @ARGV;
my $pdb_name = shift @ARGV;
undef @filenames;
$_ = $pdb_name;
s/\.pdb//;
@filenames = $_.".simplenet";
my $outindex_temp = $simple_net_index."/".$outindex_temp_tag."_net/";
if(-e $outindex_temp){
} else {
	mkdir ($outindex_temp);
}
foreach  (@filenames) {
	$inpath_temp = $simple_net_index."/".$_;
	&simplenet2ynet($inpath_temp,$outindex_temp,$_);
}
sub simplenet2ynet {
	my $filepath = $_[0];
	open(FILES,"$filepath");
	@files = <FILES>;
	close(FILES);
	$outindex = $_[1];
	$filename = $_[2];
	$_ = $filename; 
	s/\_\.simplenet//;
	$filename = $_;
	for ( my $i=0;$i<@files;$i++ ) {
		if ( 1 ) {
			$tyr_row = $i;
			@tyr_chain_site_ = split/\s+/,$files[$i];
			$tyr_chain_site = $tyr_chain_site_[0];
			@temp_outpath = split/\_/,$tyr_chain_site;
			$_ = $filename;
			s/\.simplenet//;
			if ( $temp_outpath[2] eq '' ) {
				next;
			}
			$outpath = $outindex.$_."_".$temp_outpath[1]."_".$temp_outpath[2]."_".$temp_outpath[0].".net";
			open(OUT,">$outpath");
			for ( my $j=0;$j<@files;$j++ ) {
				@temp_j_ = split/\s+/,$files[$j];
				if ( $temp_j_[0] eq $tyr_chain_site ) {
					@temp_j = split/\s+/,$files[$j];
					syswrite OUT, $temp_j[0]." ".$temp_j[1]." ".$temp_j[2]."\n";
				}
				if ( $temp_j_[2] eq $tyr_chain_site ) {
					@temp_j = split/\s+/,$files[$j];
					syswrite OUT, $temp_j[2]." ".$temp_j[1]." ".$temp_j[0]."\n";
				}
			}
			close(OUT);
			open(FILES1,"$outpath");
			@files1 = <FILES1>;
			close(FILES1);
			open(OUT2,">>$outpath");
			for ( my $ii=0;$ii<@files1;$ii++ ) {
				@temp_row1 = split/\s+/,$files1[$ii];
				$tyr_chain_site1 = $temp_row1[0];
				chomp($temp_row1[2]);
				$aa_chain_site1 = $temp_row1[2];
				for ( my $jj=0;$jj<@files;$jj++ ) {
					chomp($files[$jj]);
					@temp_jj_ = split/\s+/,$files[$jj];
					if ( ($temp_jj_[0] eq $aa_chain_site1)&&(!($temp_jj_[2] eq $tyr_chain_site1)) ) {
						chomp($files[$jj]);
						@temp_jj = split/\s+/,$files[$jj];
						syswrite OUT2, $temp_jj[0]." ".$temp_jj[1]." ".$temp_jj[2]."\n";
					}
					if ( ($temp_jj_[2] eq $aa_chain_site1)&&(!($temp_jj_[0] eq $tyr_chain_site1)) ) {
						chomp($files[$jj]);
						@temp_jj = split/\s+/,$files[$jj];
						syswrite OUT2, $temp_jj[2]." ".$temp_jj[1]." ".$temp_jj[0]."\n";
					}
				}
			}
			close(OUT);
			open(FILES2,"$outpath");
			@files2 = <FILES2>;
			close(FILES2);
			open(OUT3,">>$outpath");
			@yno2_ = split/\s+/,$files2[0];
			$yno2 = $yno2_[0];
			for ( my $iii=0;$iii<@files2;$iii++ ) {
				@temp_row2 = split/\s+/,$files2[$iii];
				if ( $temp_row2[0] ne $yno2 ) {
					$aa1_chain_site = $temp_row2[0];
					$aa2_chain_site = $temp_row2[2];
					my $tag = 0;
					for ( my $jjj=0;$jjj<@files2;$jjj++ ) {
						if ( ($files2[$jjj] =~ /$yno2/)&&($files2[$jjj] =~ /$aa2_chain_site/) ) {
							$tag++;
						}
					}
					for ( $jjj=0;$jjj<@files;$jjj++ ) {
						if ( ($files[$jjj] =~ /^$aa2_chain_site /)&&((!($files[$jjj] =~ /$yno2/))&&(!($files[$jjj] =~ /$aa1_chain_site/))&&($tag<1)) ) {
							@temp_jjj = split/\s+/,$files[$jjj];
							syswrite OUT3, $temp_jjj[0]." ".$temp_jjj[1]." ".$temp_jjj[2]."\n";
						}
						if ( ($files[$jjj] =~ / $aa2_chain_site$/)&&((!($files[$jjj] =~ /$yno2/))&&(!($files[$jjj] =~ /$aa1_chain_site/))&&($tag<1)) ) {
							@temp_jjj = split/\s+/,$files[$jjj];
							syswrite OUT3, $temp_jjj[2]." ".$temp_jjj[1]." ".$temp_jjj[0]."\n";
						}
					}
				}
			}
		}
	}
}