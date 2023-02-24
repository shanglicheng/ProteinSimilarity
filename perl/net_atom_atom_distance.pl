#!/usr/bin/perl
my $radius = 6;
my $pdb_index = shift @ARGV;
$pdb_index = $pdb_index."/";
my $pdb_file = shift @ARGV;
$_ = $pdb_file;
s/\.pdb//;
&nearest_n_residue("$pdb_index",$_,"$pdb_index",$_,$nsite);
sub nearest_n_residue {# 0:$pdbfile_index; 1:$pdb_name; 2: the number of the nearest residues 3: the out_file index of the list of the nearest residues;
	my $file_index = $_[0];
	my $file_name = $_[1];
	my $out_file_index = $_[2];
	my $protein_name_temp = $_[3];
	my $pdb_filepath = $file_index . $file_name . ".pdb";
	my $out_filepath_samechain = $out_file_index . $protein_name_temp."_".$protein_site_temp.".atom_atom";
	open (TYR_AMINO_NEAREST_FILE_SAMECHAIN,">$out_filepath_samechain");
	open (PDB_FILE,"<$pdb_filepath");
	@pdb_file = <PDB_FILE>;
	close(PDB_FILE);
	$row_num = @pdb_file;
	for ( my $i=0;$i<$row_num;$i++ ) {
		if ( $pdb_file[$i] =~ /^MODEL        1/ ) {
			$stop_row_num_temp = $i;
			last;
		} else {
			$stop_row_num_temp = 0;
		}
	}
	if ( $stop_row_num_temp != 0 ) {
		$stop_row_num = $stop_row_num_temp;
		until ( $pdb_file[$stop_row_num] =~ /^ENDMDL/ ) {
			$stop_row_num++;
		}
		$row_num = $stop_row_num;
	}#	print $row_num."\n";
	my $chain_num = 0;
	$start_chain_row[$chain_num] = 0;
	for ( my $i=0;$i<$row_num;$i++ ) {
		if ( $pdb_file[$i] =~ /^TER/ ) {
			$end_chain_row[$chain_num] = $i-1;
			$chain_num++;
			$start_chain_row[$chain_num] = $i+1;
		}
	}
	for ( my $i=$start_chain_row[0];$i<$start_chain_row[$chain_num];$i++ ) {
		if ( ($pdb_file[$i] =~ /^ATOM/)||($pdb_file[$i] =~/^HETATM/) ) {
			$row_temp_i = $pdb_file[$i];
			$residue_i = substr($row_temp_i,17,3);
			$chain_i = substr($row_temp_i,21,1);
			$resseq_i = substr($row_temp_i,22,4);
			$_ = $resseq_i;
			s/^\s+//;
			$resseq_i = $_;
			$atom_i = substr($row_temp_i,12,4);
			$_ = $atom_i;
			s/^\s+//;
			s/\s+$//;
			$atom_i = $_;
			$serial_i = substr($row_temp_i,6,5);
			$_ = $serial_i;
			s/^\s+//;
			s/\s+$//;
			$serial_i = $_;
			$x_i = substr($row_temp_i,30,8);
			$y_i = substr($row_temp_i,38,8);
			$z_i = substr($row_temp_i,46,8);
			$element_i = substr($row_temp_i,76,2);
			$residue_site_i = $residue_i.$resseq_i;
			for ( my $j=$i+1;$j<$start_chain_row[$chain_num];$j++) {
				if ( ($pdb_file[$j] =~ /^ATOM/)||($pdb_file[$j] =~/^HETATM/) ) {
					$row_temp_j = $pdb_file[$j];
					$residue_j = substr($row_temp_j,17,3);
					$chain_j = substr($row_temp_j,21,1);
					$resseq_j = substr($row_temp_j,22,4);
					$_ = $resseq_j;
					s/^\s+//;
					$resseq_j = $_;
					$atom_j = substr($row_temp_j,12,4);
					$_ = $atom_j;
					s/^\s+//;
					s/\s+$//;
					$atom_j = $_;
					$serial_j = substr($row_temp_j,6,5);
					$_ = $serial_j;
					s/^\s+//;
					s/\s+$//;
					$serial_j = $_;
					$x_j = substr($row_temp_j,30,8);
					$y_j = substr($row_temp_j,38,8);
					$z_j = substr($row_temp_j,46,8);
					$element_j = substr($row_temp_j,76,2);
					$residue_site_j = $residue_j.$resseq_j;
					if ( (($element_i =~ /c/i)||($element_i =~ /n/i)||($element_i =~ /o/i)||($element_i =~ /s/i))&&
						(($element_j =~ /c/i)||($element_j =~ /n/i)||($element_j =~ /o/i)||($element_j =~ /s/i))&&($residue_site_i ne $residue_site_j) ) {
						$dis = sqrt(($x_i-$x_j)*($x_i-$x_j)+($y_i-$y_j)*($y_i-$y_j)+($z_i-$z_j)*($z_i-$z_j));
						$dis_temp = sprintf("%.3f", $dis);
						if ( $dis_temp < $radius ) {
							syswrite TYR_AMINO_NEAREST_FILE_SAMECHAIN, $residue_i."_".$resseq_i."_".$chain_i."_".$atom_i."_".$serial_i." ";
							syswrite TYR_AMINO_NEAREST_FILE_SAMECHAIN, $residue_j."_".$resseq_j."_".$chain_j."_".$atom_j."_".$serial_j." ".$dis_temp."\n";
						}
					}
				}
			}
		}
	}
}