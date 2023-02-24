#!/usr/bin/perl
my $net_index_temp = shift @ARGV;
my $net_index_tag = shift @ARGV;
my $pdb_name = shift @ARGV;
$_ = $pdb_name;
s/\.pdb//;
my $netindex = $net_index_temp."/".$net_index_tag."_net/";
my $out_index = $net_index_temp."/".$net_index_tag."_linmatrix_";
my $aapath = $net_index_temp."\\"."aafactor";
undef @aa;
open(AAA,"$aapath");
@aa = <AAA>;
close(AAA);
undef @filenames;
	opendir(FILENAME,"$netindex");
	@filenames = readdir FILENAME;
	close(FILENAME);
	foreach  (@filenames) {
		if ( !($_ =~ /^\./) ) {
			$path = $netindex.$_;
			s/\.net//;
			$pdb_chain_site = $_;
			open(NET,"$path");
			@nets = <NET>;
			close(NET);
			undef @node;
			undef $node_num;
			@y_ = split/\s+/,$nets[0];
			$y = $y_[0];
			$node[0] = $y;
			$node_num = 1;
			for ( my $i=0;$i<@nets;$i++ ) {
				@temp_node = split/\s+/,$nets[$i];
				if ( $temp_node[0] eq $y ) {
					$node[$node_num] = $temp_node[2];
					$node_num++;
				}
			}
			$node_num_1 = $node_num;
			$node_num_temp = 0;
			for ( my $i=1;$i<$node_num_1;$i++ ) {
				for ( my $j=0;$j<@nets;$j++ ) {
					@temp_2_j = split/\s+/,$nets[$j];
					if ( $node[$i] eq $temp_2_j[0] ) {
						$node[$node_num] = $temp_2_j[2];
						$node_num++;
					}
				}
			}
			undef @in_;
			undef @out_;
			for ( my $i=0;$i<$node_num;$i++ ) {
				 push (@in_,$node[$i]);
			}
			undef %saw; 
			@out_ = grep(!$saw{$_}++,@in_);
			undef @node;
			@node = @out_;
			$node_num = @node;
			for ( my $i=0;$i<$node_num;$i++ ) {
				for ( my $j=0;$j<$node_num;$j++ ) {
					$node_matrix[$i][$j] = 0;
				}
			}
			for ( my $i=0;$i<$node_num;$i++ ) {
				for ( my $j=0;$j<$node_num;$j++ ) {
					for ( my $k=0;$k<@nets;$k++ ) {
						@temp_row = split/\s+/,$nets[$k];
						if ( ($temp_row[0] eq $node[$i])&&($temp_row[2] eq $node[$j]) ) {
							$node_matrix[$i][$j] = 1;
							$node_matrix[$j][$i] = 1;
						}
					}
				}
			}
			undef @net_temp_temp;
			@net_temp_temp = @nets;
			for ( my $i=0;$i<$node_num;$i++ ) {
				for ( my $j=0;$j<@nets;$j++ ) {
					@temp_ = split/\s+/,$nets[$j];
					if ( @temp_[2] eq $node[$i] ) {
						$temp = $temp_[2]." ".$temp_[1]." ".$temp_[0]."\n";
						push (@net_temp_temp,$temp);
					}
				}
			}
			@in__ = @net_temp_temp;
			undef %saw;
			@out__ = grep(!$saw{$_}++,@in__);
			@nets_temp = @out__;
			for ( my $i=0;$i<$node_num;$i++ ) {
				$node_degree[$i] = 0;
				for ( my $k=0;$k<@nets_temp;$k++ ) {
					@temp_node_d = split/\s+/,$nets_temp[$k];
					if ( $temp_node_d[0] eq $node[$i] ) {
						$node_degree[$i]++;
					}
				}
			}
for ( my $aa_num=0;$aa_num<@aa;$aa_num++ ) {
	my $outindex = $out_index.$aa_num."\\";
	if(-e $outindex){
	} else {
		mkdir ($outindex);
	}
	$outpath = $outindex.$pdb_chain_site.".linjie.matrix";
			open(OUT,">$outpath");
			for ( my $i=0;$i<$node_num;$i++ ) {
				@temp_aa = split/\_/,$node[$i];
				$node_aa_temp = $temp_aa[0];
				@pro = split/\s+/,$aa[$aa_num];
				$res = &residue_name($node_aa_temp);
				$col_i = &residue_property($res);
				$node_aa[$i] = $pro[$col_i];
			}
			for ( my $i=0;$i<$node_num;$i++ ) {
				$node_degree_temp = sqrt($node_degree[$i]);
				$node_degree_temp = sprintf("%.3f",$node_degree_temp);
				$node_aa_temp = sqrt($node_aa[$i]);
				$node_aa_temp = sprintf("%.3f",$node_aa_temp);
				syswrite OUT, $node_degree_temp." ".$node_aa_temp." ";
				for ( my $j=0;$j<$node_num;$j++ ) {
					if ( $j < ($node_num-1) ) {
						syswrite OUT, $node_matrix[$i][$j]." ";
					} else {
						syswrite OUT, $node_matrix[$i][$j];
					}
				}
				syswrite OUT, "\n";
			}
		}
	}
}
sub residue_name {
	my $sub_residue = $_[0];
	my $residue_pro = ' ';
	if ($sub_residue =~ /ala/i) {
		$residue_pro = A;	}
	if ($sub_residue =~ /sbd/i) {
		$residue_pro = A;	}
	if ($sub_residue =~ /maa/i) {
		$residue_pro = A;	}
	if ($sub_residue =~ /dal/i) {
		$residue_pro = A;	}
	if ($sub_residue =~ /clb/i) {
		$residue_pro = A;	}
	if ($sub_residue =~ /aba/i) {
		$residue_pro = A;	}
	if ($sub_residue =~ /sec/i) {
		$residue_pro = A;	}
	if ($sub_residue =~ /dab/i) {
		$residue_pro = A;	}
	if ($sub_residue =~ /cld/i) {
		$residue_pro = A;	}
	if ($sub_residue =~ /dam/i) {
		$residue_pro = A;	}
	if ($sub_residue =~ /arg/i) {
		$residue_pro = R;	}
	if ($sub_residue =~ /asn/i) {
		$residue_pro = N;	}
	if ($sub_residue =~ /asx/i) {
		$residue_pro = N;	}
	if ($sub_residue =~ /asp/i) {
		$residue_pro = D;	}
	if ($sub_residue =~ /ias/i) {
		$residue_pro = D;	}
	if ($sub_residue =~ /cys/i) {
		$residue_pro = C;	}
	if ($sub_residue =~ /ocs/i) {
		$residue_pro = C;	}
	if ($sub_residue =~ /cso/i) {
		$residue_pro = C;	}
	if ($sub_residue =~ /sch/i) {
		$residue_pro = C;	}
	if ($sub_residue =~ /csx/i) {
		$residue_pro = C;	}
	if ($sub_residue =~ /ccs/i) {
		$residue_pro = C;	}
	if ($sub_residue =~ /cme/i) {
		$residue_pro = C;	}
	if ($sub_residue =~ /cas/i) {
		$residue_pro = C;	}
	if ($sub_residue =~ /cmt/i) {
		$residue_pro = C;	}
	if ($sub_residue =~ /scy/i) {
		$residue_pro = C;	}
	if ($sub_residue =~ /gln/i) {
		$residue_pro = Q;	}
	if ($sub_residue =~ /glu/i) {
		$residue_pro = E;	}
	if ($sub_residue =~ /dgl/i) {
		$residue_pro = E;	}
	if ($sub_residue =~ /gly/i) {
		$residue_pro = G;	}
	if ($sub_residue =~ /his/i) {
		$residue_pro = H;	}
	if ($sub_residue =~ /hic/i) {
		$residue_pro = H;	}
	if ($sub_residue =~ /ile/i) {
		$residue_pro = I;	}
	if ($sub_residue =~ /leu/i) {
		$residue_pro = L;	}
	if ($sub_residue =~ /lef/i) {
		$residue_pro = L;	}
	if ($sub_residue =~ /lys/i) {
		$residue_pro = K;	}
	if ($sub_residue =~ /xx1/i) {
		$residue_pro = K;	}
	if ($sub_residue =~ /dmo/i) {
		$residue_pro = K;	}
	if ($sub_residue =~ /mly/i) {
		$residue_pro = K;	}
	if ($sub_residue =~ /llp/i) {
		$residue_pro = K;	}
	if ($sub_residue =~ /met/i) {
		$residue_pro = M;	}
	if ($sub_residue =~ /mse/i) {
		$residue_pro = M;	}
	if ($sub_residue =~ /phe/i) {
		$residue_pro = F;	}
	if ($sub_residue =~ /mea/i) {
		$residue_pro = F;	}
	if ($sub_residue =~ /pro/i) {
		$residue_pro = P;	}
	if ($sub_residue =~ /ser/i) {
		$residue_pro = S;	}
	if ($sub_residue =~ /dsn/i) {
		$residue_pro = S;	}
	if ($sub_residue =~ /sbl/i) {
		$residue_pro = S;	}
	if ($sub_residue =~ /sep/i) {
		$residue_pro = S;	}
	if ($sub_residue =~ /thr/i) {
		$residue_pro = T;	}
	if ($sub_residue =~ /tpo/i) {
		$residue_pro = T;	}
	if ($sub_residue =~ /trp/i) {
		$residue_pro = W;	}
	if ($sub_residue =~ /ftr/i) {
		$residue_pro = W;	}
	if ($sub_residue =~ /tyr/i) {
		$residue_pro = Y;	}
	if ($sub_residue =~ /ynm/i) {
		$residue_pro = Y;	}
	if ($sub_residue =~ /tys/i) {
		$residue_pro = Y;	}
	if ($sub_residue =~ /ptr/i) {
		$residue_pro = Y;	}
	if ($sub_residue =~ /niy/i) {
		$residue_pro = Y;	}
	if ($sub_residue =~ /yof/i) {
		$residue_pro = Y;	}
	if ($sub_residue =~ /val/i) {
		$residue_pro = V;	}
	$residue_pro;
}
sub residue_property {
	my $sub_residue = $_[0];
	$residue_pro = 21;
	if ($sub_residue =~ /A/i) {
		$residue_pro = 0;
	}
	if ($sub_residue =~ /R/i) {
		$residue_pro = 1;
	}
	if ($sub_residue =~ /N/i) {
		$residue_pro = 2;
	}
	if ($sub_residue =~ /D/i) {
		$residue_pro = 3;	
	}
	if ($sub_residue =~ /C/i) {
		$residue_pro = 4;
	}
	if ($sub_residue =~ /Q/i) {		
		$residue_pro = 5;	
	}
	if ($sub_residue =~ /E/i) {		
		$residue_pro = 6;	
	}
	if ($sub_residue =~ /G/i) {
		$residue_pro = 7;	
	}
	if ($sub_residue =~ /H/i) {
		$residue_pro = 8;	
	}
	if ($sub_residue =~ /I/i) {
		$residue_pro = 9;	
	}
	if ($sub_residue =~ /L/i) {
		$residue_pro = 10;	
	}
	if ($sub_residue =~ /K/i) {
		$residue_pro = 11;	
	}
	if ($sub_residue =~ /M/i) {
		$residue_pro = 12;	
	}
	if ($sub_residue =~ /F/i) {
		$residue_pro = 13;	
	}
	if ($sub_residue =~ /P/i) {	
		$residue_pro = 14;	
	}
	if ($sub_residue =~ /S/i) {
		$residue_pro = 15;	
	}
	if ($sub_residue =~ /T/i) {
		$residue_pro = 16;	
	}
	if ($sub_residue =~ /W/i) {	
		$residue_pro = 17;
	}
	if ($sub_residue =~ /Y/i) {	
		$residue_pro = 18;	
	}
	if ($sub_residue =~ /V/i) {
		$residue_pro = 19;	
	}
	$residue_pro;
}