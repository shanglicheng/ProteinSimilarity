#!/usr/bin/perl
my $pdb_index = shift @ARGV;
$model_txt_path = $pdb_index."/model.txt";
$targ_txt_path = $pdb_index."/targ.txt";
$test_txt_path = $pdb_index."/test.txt";
my $outpath = $pdb_index."/result.out";
print "################################################\n";
print "\tWelcome to use Protein Analysis Software\n";
print "\tSoftware working directory is $pdb_index\n";
print "################################################\n";
print "Running........";
open(OUT,">$outpath");
open(FILE,"$model_txt_path");
@pdb_models = <FILE>;
close(FILE);
open(FILE,"$targ_txt_path");
@pdb_targs = <FILE>;
close(FILE);
open(FILE,"$test_txt_path");
@pdb_tests = <FILE>;
close(FILE);
undef @pdb_;
foreach  (@pdb_models) {
	chomp($_);
	@temp = split/\s+/,$_;
	push @pdb_,$temp[0];
}
undef @in_;
undef @out_;
for ( my $i=0;$i<@pdb_;$i++ ) {
	 push (@in_,$pdb_[$i]);
}
undef %saw; 
@out_ = grep(!$saw{$_}++,@in_);
@pdb = @out_;
foreach  (@pdb) {
	$pdb_model = $_.".pdb";
	system("net_atom_atom_distance.pl $pdb_index $pdb_model");
	system("net_atom_by_vanderwaal.pl $pdb_index $pdb_model");
	system("net_residue_by_vanderwaal.pl $pdb_index $pdb_model");
	system("net_res_simple_band.pl $pdb_index $pdb_model");
	system("all_residue_net_2_Y_net.pl $pdb_index model $pdb_model");
	system("resnet2linjie_matrix_2level_1distance.pl $pdb_index model $pdb_model");
	system("cal_eig.pl $pdb_index model $pdb_model");
}
foreach  (@pdb_tests) {
	chomp($_);
	@temp = split/\s+/,$_;
	$pdb_targ = $temp[0].".pdb";
	system("net_atom_atom_distance.pl $pdb_index $pdb_targ");
	system("net_atom_by_vanderwaal.pl $pdb_index $pdb_targ");
	system("net_residue_by_vanderwaal.pl $pdb_index $pdb_targ");
	system("net_res_simple_band.pl $pdb_index $pdb_targ");
	system("all_residue_net_2_Y_net.pl $pdb_index test $pdb_targ");
	system("resnet2linjie_matrix_2level_1distance.pl $pdb_index test $pdb_targ");
	system("cal_eig.pl $pdb_index test $pdb_targ");
}
foreach  (@pdb_targs) {
	chomp($_);
	@temp = split/\s+/,$_;
	$pdb_targ = $temp[0].".pdb";
	system("net_atom_atom_distance.pl $pdb_index $pdb_targ");
	system("net_atom_by_vanderwaal.pl $pdb_index $pdb_targ");
	system("net_residue_by_vanderwaal.pl $pdb_index $pdb_targ");
	system("net_res_simple_band.pl $pdb_index $pdb_targ");
	system("all_residue_net_2_Y_net.pl $pdb_index targ $pdb_targ");
	system("resnet2linjie_matrix_2level_1distance.pl $pdb_index targ $pdb_targ");
	system("cal_eig.pl $pdb_index targ $pdb_targ");
}
$cindex_model = $pdb_index."/model_character_index/";
$cindex_test = $pdb_index."/test_character_index/";
$cindex_targ = $pdb_index."/targ_character_index/";
opendir(FILENAME,"$cindex_model");
@model_files = readdir FILENAME;
close(FILENAME);
$model_num = 0;
undef @model_array;
foreach  (@pdb_models) {
	chomp($_);
	@temp = split/\s+/,$_;
	$model_name = $temp[0];
	$pdb_model_chain = $temp[1];
	$model_site_ = $temp[2];
	for ( my $i=0;$i<@model_files;$i++ ) {
		if ( !($model_files[$i] =~ /^\./) ) {
			$model_file_name = $model_files[$i];
			$_ = $model_file_name;
			s/\.cindex//;
			@temp_ = split/\_/,$_;
			if ( $temp_[3] == $model_site_ ) {
				$filename = $model_file_name;
			}
		}
	}
	$model_path = $cindex_model.$filename;
	$model_protein_chain_site = $model_name."_".$pdb_model_chain."_".$model_site_;
	open(FILE,"$model_path");
	@a_model = <FILE>;
	close(FILE);
	$model_array[$model_num] = "";
	for ( my $i=0;$i<@a_model;$i++ ) {
		$a_model_temp = $a_model[$i];
		chomp($a_model_temp);
		$model_array[$model_num] = $model_array[$model_num]." ".$a_model_temp;
	}
	$_ = $model_array[$model_num];
	s/^\s+//;
	s/\s+$//;
	$model_array[$model_num] = $_;
	$model_array_tag[$model_num] = $model_protein_chain_site;
	$model_num++;
}
opendir(FILENAME,"$cindex_test");
@test_files = readdir FILENAME;
close(FILENAME);
$test_num = 0;
undef @test_array;
foreach  (@pdb_tests) {
	chomp($_);
	@temp = split/\s+/,$_;
	$test_name = $temp[0];
	$pdb_test_chain = $temp[1];
	$test_site_ = $temp[2];
	for ( my $i=0;$i<@test_files;$i++ ) {
		if ( !($test_files[$i] =~ /^\./) ) {
			$test_file_name = $test_files[$i];
			$_ = $test_file_name;
			s/\.cindex//;
			@temp_ = split/\_/,$_;
			if ( $temp_[3] == $test_site_ ) {
				$filename = $test_file_name;
			}
		}
	}
	$test_path = $cindex_test.$filename;
	$test_protein_chain_site = $test_name."_".$pdb_test_chain."_".$test_site_;
	open(FILE,"$test_path");
	@a_test = <FILE>;
	close(FILE);
	$test_array[$test_num] = "";
	for ( my $i=0;$i<@a_test;$i++ ) {
		$a_test_temp = $a_test[$i];
		chomp($a_test_temp);
		$test_array[$test_num] = $test_array[$test_num]." ".$a_test_temp;
	}
	$_ = $test_array[$test_num];
	s/^\s+//;
	s/\s+$//;
	$test_array[$test_num] = $_;
	$test_array_tag[$test_num] = $test_protein_chain_site;
	$test_num++;
}
undef @acu2;
undef @sen2;
undef @mcc2;
for ( my $select1=0;$select1<5 ;$select1++) {
for ( my $i=0;$i<5;$i++ ) {
	if ( $i != $select1 ) {
	$tp = 0;
	$fp = 0;
	$tn = 0;
	$fn = 0;
	for ( my $ii=0;$ii<$model_num;$ii++ ) {
		undef @sim_model;
		$sim_model_num = 0;
		for ( my $iii=0;$iii<$model_num ;$iii++) {
			if ( $iii != $ii ) {
				@temp_model_array_ = split/\s+/,$model_array[$ii];
				@temp_test_array_ = split/\s+/,$model_array[$iii];
				$temp_model_array = $temp_model_array_[$select1]." ".$temp_model_array_[$i];
				$temp_test_array = $temp_test_array_[$select1]." ".$temp_test_array_[$i];
				$sim_model[$sim_model_num] = &array_sim($temp_model_array,$temp_test_array);
				$sim_model_num++;
			}
		}
		$sim_test_num = 0;
		for ( my $iii=0;$iii<$test_num ;$iii++) {
				@temp_model_array_ = split/\s+/,$model_array[$ii];
				@temp_test_array_ = split/\s+/,$test_array[$iii];
				$temp_model_array = $temp_model_array_[$select1]." ".$temp_model_array_[$i];
				$temp_test_array = $temp_test_array_[$select1]." ".$temp_test_array_[$i];
				$sim_test[$sim_test_num] = &array_sim($temp_model_array,$temp_test_array);
				$sim_test_num++;
		}
		use List::Util qw/max min/;
		$max_model = max(@sim_model);
		$max_test = max(@sim_test);
		if ( $max_model >= $max_test ) {
			$tp++;
		} else {
			$fn++;
		}
	}
	for ( my $ii=0;$ii<$test_num;$ii++ ) {
		undef @sim_model;
		$sim_model_num = 0;
		for ( my $iii=0;$iii<$model_num ;$iii++) {
				@temp_model_array_ = split/\s+/,$test_array[$ii];
				@temp_test_array_ = split/\s+/,$model_array[$iii];
				$temp_model_array = $temp_model_array_[$select1]." ".$temp_model_array_[$i];
				$temp_test_array = $temp_test_array_[$select1]." ".$temp_test_array_[$i];
				$sim_model[$sim_model_num] = &array_sim($temp_model_array,$temp_test_array);
				$sim_model_num++;
		}
		$sim_test_num = 0;
		for ( my $iii=0;$iii<$test_num ;$iii++) {
			if ( $iii != $ii ) {
				@temp_model_array_ = split/\s+/,$test_array[$ii];
				@temp_test_array_ = split/\s+/,$test_array[$iii];
				$temp_model_array = $temp_model_array_[$select1]." ".$temp_model_array_[$i];
				$temp_test_array = $temp_test_array_[$select1]." ".$temp_test_array_[$i];
				$sim_test[$sim_test_num] = &array_sim($temp_model_array,$temp_test_array);
				$sim_test_num++;
			}
		}
		use List::Util qw/max min/;
		$max_model = max(@sim_model);
		$max_test = max(@sim_test);
		if ( $max_model >= $max_test ) {
			$fp++;
		} else {
			$tn++;
		}
	}
	$acu2[$select1][$i] = ($tp+$tn)/($tp+$fp+$tn+$fn);
	$sen2[$select1][$i] = $tp/($tp+$fp);
	$mcc2[$select1][$i] = abs(($tp*$tn-$fp*$fn)/(sqrt(($tp+$fp)*($tp+$fn)*($tn+$fp)*($tn+$fn))));
	} else {
		$mcc2[$select1][$i] = 0;
	}
}
}
$max = 0;
for ( my $se=0;$se<5;$se++ ) {
	for (my $see=0;$see<5;$see++ ) {
		if ( $max < $mcc2[$se][$see] ) {
			$max = $mcc2[$se][$see];
			$max_col1 = $se;
			$max_col2 = $see;
			$max_auc2 = $acu2[$se][$see];
			$max_sen2 = $sen2[$se][$see];
		}
	}
}
$select1 = $max_col1;
$select2 = $max_col2;
$auc2_ = $max_auc2;
$sen2_ = $max_sen2;
undef @acu3;
undef @sen3;
undef @mcc3;
for ( my $i=0;$i<5;$i++ ) {
	if ( ($i != $select1)&&($i != $select2) ) {
	$tp = 0;
	$fp = 0;
	$tn = 0;
	$fn = 0;
	for ( my $ii=0;$ii<$model_num;$ii++ ) {
		undef @sim_model;
		$sim_model_num = 0;
		for ( my $iii=0;$iii<$model_num ;$iii++) {
			if ( $iii != $ii ) {
				@temp_model_array_ = split/\s+/,$model_array[$ii];
				@temp_test_array_ = split/\s+/,$model_array[$iii];
				$temp_model_array = $temp_model_array_[$select1]." ".$temp_model_array_[$select2]." ".$temp_model_array_[$i];
				$temp_test_array = $temp_test_array_[$select1]." ".$temp_test_array_[$select2]." ".$temp_test_array_[$i];
				$sim_model[$sim_model_num] = &array_sim($temp_model_array,$temp_test_array);
				$sim_model_num++;
			}
		}
		$sim_test_num = 0;
		for ( my $iii=0;$iii<$test_num ;$iii++) {
				@temp_model_array_ = split/\s+/,$model_array[$ii];
				@temp_test_array_ = split/\s+/,$test_array[$iii];
				$temp_model_array = $temp_model_array_[$select1]." ".$temp_model_array_[$select2]." ".$temp_model_array_[$i];
				$temp_test_array = $temp_test_array_[$select1]." ".$temp_test_array_[$select2]." ".$temp_test_array_[$i];
				$sim_test[$sim_test_num] = &array_sim($temp_model_array,$temp_test_array);
				$sim_test_num++;
		}
		use List::Util qw/max min/;
		$max_model = max(@sim_model);
		$max_test = max(@sim_test);
		if ( $max_model >= $max_test ) {
			$tp++;
		} else {
			$fp++;
		}
	}
	for ( my $ii=0;$ii<$test_num;$ii++ ) {
		undef @sim_model;
		$sim_model_num = 0;
		for ( my $iii=0;$iii<$model_num ;$iii++) {
			
				@temp_model_array_ = split/\s+/,$test_array[$ii];
				@temp_test_array_ = split/\s+/,$model_array[$iii];
				$temp_model_array = $temp_model_array_[$select1]." ".$temp_model_array_[$select2]." ".$temp_model_array_[$i];
				$temp_test_array = $temp_test_array_[$select1]." ".$temp_test_array_[$select2]." ".$temp_test_array_[$i];
				$sim_model[$sim_model_num] = &array_sim($temp_model_array,$temp_test_array);
				$sim_model_num++;
			
		}
		$sim_test_num = 0;
		for ( my $iii=0;$iii<$test_num ;$iii++) {
			if ( $iii != $ii ) {
				@temp_model_array_ = split/\s+/,$test_array[$ii];
				@temp_test_array_ = split/\s+/,$test_array[$iii];
				$temp_model_array = $temp_model_array_[$select1]." ".$temp_model_array_[$select2]." ".$temp_model_array_[$i];
				$temp_test_array = $temp_test_array_[$select1]." ".$temp_test_array_[$select2]." ".$temp_test_array_[$i];
				$sim_test[$sim_test_num] = &array_sim($temp_model_array,$temp_test_array);
				$sim_test_num++;
			}
		}
		use List::Util qw/max min/;
		$max_model = max(@sim_model);
		$max_test = max(@sim_test);
		if ( $max_model >= $max_test ) {
			$fn++;
		} else {
			$tn++;
		}
	}
	
	$acu3[$i] = ($tp+$tn)/($tp+$fp+$tn+$fn);
	$sen3[$i] = $tp/($tp+$fp);
	$mcc3[$i] = abs(($tp*$tn-$fp*$fn)/(sqrt(($tp+$fp)*($tp+$fn)*($tn+$fp)*($tn+$fn))));
	} else {
		$mcc3[$i] = 0;
	}
}
use List::Util qw/max min/;
$max_mcc3 = max(@mcc3);
for ( my $k=0;$k<@mcc3;$k++ ) {
	if ( $max_mcc3 == $mcc3[$k] ) {
		$select3 = $k;
		$mcc3_ = $mcc3[$k];
		$acu3_ = $acu3[$k];
		$sen3_ = $sen3[$k];
	}
}
undef @acu4;
undef @sen4;
undef @mcc4;
for ( my $i=0;$i<5;$i++ ) {
	if ( ($i != $select1)&&($i != $select2)&&($i != $select3) ) {
	$tp = 0;
	$fp = 0;
	$tn = 0;
	$fn = 0;
	for ( my $ii=0;$ii<$model_num;$ii++ ) {
		undef @sim_model;
		$sim_model_num = 0;
		for ( my $iii=0;$iii<$model_num ;$iii++) {
			if ( $iii != $ii ) {
				@temp_model_array_ = split/\s+/,$model_array[$ii];
				@temp_test_array_ = split/\s+/,$model_array[$iii];
				$temp_model_array = $temp_model_array_[$select1]." ".$temp_model_array_[$select2]." ".$temp_model_array_[$select3]." ".$temp_model_array_[$i];
				$temp_test_array = $temp_test_array_[$select1]." ".$temp_test_array_[$select2]." ".$temp_test_array_[$select3]." ".$temp_test_array_[$i];
				$sim_model[$sim_model_num] = &array_sim($temp_model_array,$temp_test_array);
				$sim_model_num++;
			}
		}
		$sim_test_num = 0;
		for ( my $iii=0;$iii<$test_num ;$iii++) {
				@temp_model_array_ = split/\s+/,$model_array[$ii];
				@temp_test_array_ = split/\s+/,$test_array[$iii];
				$temp_model_array = $temp_model_array_[$select1]." ".$temp_model_array_[$select2]." ".$temp_model_array_[$select3]." ".$temp_model_array_[$i];
				$temp_test_array = $temp_test_array_[$select1]." ".$temp_test_array_[$select2]." ".$temp_test_array_[$select3]." ".$temp_test_array_[$i];
				$sim_test[$sim_test_num] = &array_sim($temp_model_array,$temp_test_array);
				$sim_test_num++;
		}
		use List::Util qw/max min/;
		$max_model = max(@sim_model);
		$max_test = max(@sim_test);
		if ( $max_model >= $max_test ) {
			$tp++;
		} else {
			$fp++;
		}
	}
	for ( my $ii=0;$ii<$test_num;$ii++ ) {
		undef @sim_model;
		$sim_model_num = 0;
		for ( my $iii=0;$iii<$model_num ;$iii++) {
				@temp_model_array_ = split/\s+/,$test_array[$ii];
				@temp_test_array_ = split/\s+/,$model_array[$iii];
				$temp_model_array = $temp_model_array_[$select1]." ".$temp_model_array_[$select2]." ".$temp_model_array_[$select3]." ".$temp_model_array_[$i];
				$temp_test_array = $temp_test_array_[$select1]." ".$temp_test_array_[$select2]." ".$temp_test_array_[$select3]." ".$temp_test_array_[$i];
				$sim_model[$sim_model_num] = &array_sim($temp_model_array,$temp_test_array);
				$sim_model_num++;
		}
		$sim_test_num = 0;
		for ( my $iii=0;$iii<$test_num ;$iii++) {
			if ( $iii != $ii ) {
				@temp_model_array_ = split/\s+/,$test_array[$ii];
				@temp_test_array_ = split/\s+/,$test_array[$iii];
				$temp_model_array = $temp_model_array_[$select1]." ".$temp_model_array_[$select2]." ".$temp_model_array_[$select3]." ".$temp_model_array_[$i];
				$temp_test_array = $temp_test_array_[$select1]." ".$temp_test_array_[$select2]." ".$temp_test_array_[$select3]." ".$temp_test_array_[$i];
				$sim_test[$sim_test_num] = &array_sim($temp_model_array,$temp_test_array);
				$sim_test_num++;
			}
		}
		use List::Util qw/max min/;
		$max_model = max(@sim_model);
		$max_test = max(@sim_test);
		if ( $max_model >= $max_test ) {
			$fn++;
		} else {
			$tn++;
		}
	}
	$acu4[$i] = ($tp+$tn)/($tp+$fp+$tn+$fn);
	$sen4[$i] = $tp/($tp+$fp);
	$mcc4[$i] = abs(($tp*$tn-$fp*$fn)/(sqrt(($tp+$fp)*($tp+$fn)*($tn+$fp)*($tn+$fn))));
	} else {
		$mcc4[$i] = 0;
	}
}
use List::Util qw/max min/;
$max_mcc4 = max(@mcc4);
for ( my $k=0;$k<@mcc4;$k++ ) {
	if ( $max_mcc4 == $mcc4[$k] ) {
		$select4 = $k;
		$mcc4_ = $mcc4[$k];
		$acu4_ = $acu4[$k];
		$sen4_ = $sen4[$k];
	}
}
undef @acu5;
undef @sen5;
undef @mcc5;
for ( my $i=0;$i<5;$i++ ) {
	if ( ($i != $select1)&&($i != $select2)&&($i != $select3)&&($i != $select4) ) {
	$i_temp = $i;
	$tp = 0;
	$fp = 0;
	$tn = 0;
	$fn = 0;
	for ( my $ii=0;$ii<$model_num;$ii++ ) {
		undef @sim_model;
		$sim_model_num = 0;
		for ( my $iii=0;$iii<$model_num ;$iii++) {
			if ( $iii != $ii ) {
				@temp_model_array_ = split/\s+/,$model_array[$ii];
				@temp_test_array_ = split/\s+/,$model_array[$iii];
				$temp_model_array = $temp_model_array_[$select1]." ".$temp_model_array_[$select2]." ".$temp_model_array_[$select3]." ".$temp_model_array_[$select4]." ".$temp_model_array_[$i];
				$temp_test_array = $temp_test_array_[$select1]." ".$temp_test_array_[$select2]." ".$temp_test_array_[$select3]." ".$temp_test_array_[$select4]." ".$temp_test_array_[$i];
				$sim_model[$sim_model_num] = &array_sim($temp_model_array,$temp_test_array);
				$sim_model_num++;
			}
		}
		$sim_test_num = 0;
		for ( my $iii=0;$iii<$test_num ;$iii++) {
				@temp_model_array_ = split/\s+/,$model_array[$ii];
				@temp_test_array_ = split/\s+/,$test_array[$iii];
				$temp_model_array = $temp_model_array_[$select1]." ".$temp_model_array_[$select2]." ".$temp_model_array_[$select3]." ".$temp_model_array_[$select4]." ".$temp_model_array_[$i];
				$temp_test_array = $temp_test_array_[$select1]." ".$temp_test_array_[$select2]." ".$temp_test_array_[$select3]." ".$temp_test_array_[$select4]." ".$temp_test_array_[$i];
				$sim_test[$sim_test_num] = &array_sim($temp_model_array,$temp_test_array);
				$sim_test_num++;
		}
		use List::Util qw/max min/;
		$max_model = max(@sim_model);
		$max_test = max(@sim_test);
		if ( $max_model >= $max_test ) {
			$tp++;
		} else {
			$fp++;
		}
	}
	for ( my $ii=0;$ii<$test_num;$ii++ ) {
		undef @sim_model;
		$sim_model_num = 0;
		for ( my $iii=0;$iii<$model_num ;$iii++) {
				@temp_model_array_ = split/\s+/,$test_array[$ii];
				@temp_test_array_ = split/\s+/,$model_array[$iii];
				$temp_model_array = $temp_model_array_[$select1]." ".$temp_model_array_[$select2]." ".$temp_model_array_[$select3]." ".$temp_model_array_[$select4]." ".$temp_model_array_[$i];
				$temp_test_array = $temp_test_array_[$select1]." ".$temp_test_array_[$select2]." ".$temp_test_array_[$select3]." ".$temp_test_array_[$select4]." ".$temp_test_array_[$i];
				$sim_model[$sim_model_num] = &array_sim($temp_model_array,$temp_test_array);
				$sim_model_num++;
		}
		$sim_test_num = 0;
		for ( my $iii=0;$iii<$test_num ;$iii++) {
			if ( $iii != $ii ) {
				@temp_model_array_ = split/\s+/,$test_array[$ii];
				@temp_test_array_ = split/\s+/,$test_array[$iii];
				$temp_model_array = $temp_model_array_[$select1]." ".$temp_model_array_[$select2]." ".$temp_model_array_[$select3]." ".$temp_model_array_[$select4]." ".$temp_model_array_[$i];
				$temp_test_array = $temp_test_array_[$select1]." ".$temp_test_array_[$select2]." ".$temp_test_array_[$select3]." ".$temp_test_array_[$select4]." ".$temp_test_array_[$i];
				$sim_test[$sim_test_num] = &array_sim($temp_model_array,$temp_test_array);
				$sim_test_num++;
			}
		}
		use List::Util qw/max min/;
		$max_model = max(@sim_model);
		$max_test = max(@sim_test);
		if ( $max_model >= $max_test ) {
			$fn++;
		} else {
			$tn++;
		}
	}
	$acu5[$i] = ($tp+$tn)/($tp+$fp+$tn+$fn);
	$sen5[$i] = $tp/($tp+$fp);
	$mcc5[$i] = abs(($tp*$tn-$fp*$fn)/(sqrt(($tp+$fp)*($tp+$fn)*($tn+$fp)*($tn+$fn))));
	} else {
		$mcc5[$i] = 0;
	}
}
use List::Util qw/max min/;
$max_mcc5 = max(@mcc5);
$select5 = $i_temp;
$mcc5_ = max(@mcc5);
$acu5_ = max(@acu5);
$sen5_ = max(@sen5);
undef @mcc_all;
push @mcc_all,$mcc1_,$mcc2_,$mcc3_,$mcc4_,$mcc5_;
use List::Util qw/max min/;
$mcc_all_max = max(@mcc_all);
if ( $mcc2_ == $mcc_all_max ) {
	$auc_ = $auc2_;
	$sen_ = $sen2_;
	$select_col = $select1." ".$select2;
} elsif ( $mcc3_ == $mcc_all_max ) {
	$auc_ = $acu3_;
	$sen_ = $sen3_;
	$select_col = $select1." ".$select2." ".$select3;
} elsif ( $mcc4_ == $mcc_all_max ) {
	$auc_ = $acu4_;
	$sen_ = $sen4_;
	$select_col = $select1." ".$select2." ".$select3." ".$select4;
} elsif ( $mcc5_ == $mcc_all_max ) {
	$auc_ = $acu5_;
	$sen_ = $sen5_;
	$select_col = $select1." ".$select2." ".$select3." ".$select4." ".$select5;
}
syswrite OUT, "select property : $select_col\n";
syswrite OUT, "accuracy : $auc_\n";
syswrite OUT, "sensitivity : $sen_\n";
opendir(FILENAME,"$cindex_targ");
@targ_files = readdir FILENAME;
close(FILENAME);
$targ_num = 0;
undef @targ_array;
foreach  (@targ_files) {
	if ( !($_ =~ /^\./) ) {
	chomp($_);
	$filename = $_;
	s/\.cindex//;
	@temp = split/\_/,$_;
	$targ_name = $temp[0];
	$pdb_targ_chain = $temp[1];
	$targ_site_ = $temp[2];
	for ( my $t=0;$t<@pdb_targs;$t++) {
	$temp_tt = $pdb_targs[$t];
	chomp($temp_tt);
	@temp_t = split/\s+/,$temp_tt;
	if ( ($temp_t[1] eq $pdb_targ_chain)&&($temp_t[2] == $targ_site_)) {
	$targ_path = $cindex_targ.$filename;
	$targ_protein_chain_site = $targ_name."_".$pdb_targ_chain."_".$targ_site_;
	open(FILE,"$targ_path");
	@a_targ = <FILE>;
	close(FILE);
	$targ_array[$targ_num] = "";
	for ( my $i=0;$i<@a_targ;$i++ ) {
		$a_targ_temp = $a_targ[$i];
		chomp($a_targ_temp);
		$targ_array[$targ_num] = $targ_array[$targ_num]." ".$a_targ_temp;
	}
	$_ = $targ_array[$targ_num];
	s/^\s+//;
	s/\s+$//;
	$targ_array[$targ_num] = $_;
	$targ_array_tag[$targ_num] = $targ_protein_chain_site;
	$targ_num++;
	}
}
}
}
undef @targ_model_test;
my $num = 0;
for ( my $i=0;$i<$model_num;$i++ ) {
	for ( my $ii=0;$ii<$targ_num;$ii++ ) {
		$targ_model_test[$num][0] = "model_".$model_array_tag[$i]."_".$targ_array_tag[$ii];
		@select_cols = split/\s+/,$select_col;
		@model_array_ = split/\s+/,$model_array[$i];
		@targ_array_ = split/\s+/,$targ_array[$ii];
		$model_array_temp = "";
		$targ_array_temp = "";
		foreach  (@select_cols) {
			$model_array_temp = $model_array_temp." ".$model_array_[$_];
			$targ_array_temp = $targ_array_temp." ".$targ_array_[$_];
		}
		$_ = $model_array_temp;
		s/^\s+//;
		$model_array_temp = $_;
		$_ = $targ_array_temp;
		s/^\s+//;
		$targ_array_temp = $_;
		$targ_model_test[$num][1] = &array_sim($model_array_temp,$targ_array_temp);
		$num++;
	}
}
for ( my $i=0;$i<$test_num;$i++ ) {
	for ( my $ii=0;$ii<$targ_num;$ii++ ) {
		$targ_model_test[$num][0] = "model_".$test_array_tag[$i]."_".$targ_array_tag[$ii];
		@select_cols = split/\s+/,$select_col;
		@test_array_ = split/\s+/,$test_array[$i];
		@targ_array_ = split/\s+/,$targ_array[$ii];
		$test_array_temp = "";
		$targ_array_temp = "";
		foreach  (@select_cols) {
			$test_array_temp = $test_array_temp." ".$test_array_[$_];
			$targ_array_temp = $targ_array_temp." ".$targ_array_[$_];
		}
		$_ = $test_array_temp;
		s/^\s+//;
		$test_array_temp = $_;
		$_ = $targ_array_temp;
		s/^\s+//;
		$targ_array_temp = $_;
		$targ_model_test[$num][1] = &array_sim($test_array_temp,$targ_array_temp);
		$num++;
	}
}
syswrite OUT, "\n";
undef @re_new;
for ( my $ss=0;$ss<$num-1;$ss++ ) {
	for ( my $sss=$ss+1;$sss<$num;$sss++ ) {
		if ( $targ_model_test[$ss][1] < $targ_model_test[$sss][1] ) {
			$temp0 = $targ_model_test[$ss][0];
			$temp1 = $targ_model_test[$ss][1];
			$targ_model_test[$ss][0] = $targ_model_test[$sss][0];
			$targ_model_test[$ss][1] = $targ_model_test[$sss][1];
			$targ_model_test[$sss][0] = $temp0;
			$targ_model_test[$sss][1] = $temp1;
		}
	}
}
for ( my $i=0;$i<$num ;$i++) {
	for ( my $j=0;$j<2 ;$j++) {
		syswrite OUT, $targ_model_test[$i][$j]." ";
	}
	syswrite OUT, "\n";
}
sub array_sim{
	$a = $_[0];
	$b = $_[1];
	@a_ = split/\s+/,$a;
	@b_ = split/\s+/,$b;
	$ab = 0;
	$aa = 0;
	$bb = 0;
	for (my $i=0;$i<@a_;$i++ ) {
		$ab = $ab + $a_[$i]*$b_[$i];
		$aa = $aa + $a_[$i]*$a_[$i];
		$bb = $bb + $b_[$i]*$b_[$i];
	}
	$aa = sqrt($aa);
	$bb = sqrt($bb);
	$sim = $ab/($aa*$bb);
}
