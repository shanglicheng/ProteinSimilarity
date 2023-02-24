#!/usr/bin/perl
my $atom_atom_index = shift @ARGV;
my $atom_vannet_index = $atom_atom_index;
my $pdb_name = shift @ARGV;
my $C = 1.70;
my $N = 1.55;
my $O = 1.52;
my $S = 1.80;
my $H = 1.12;
$atom_dis_index = $atom_atom_index."/";
$out_index = $atom_vannet_index;
$_ = $pdb_name;
s/\.pdb//;
@files = $_."_.atom_atom";
foreach  (@files) {
	if ( !($_ =~ /^\./) ) {
		$filepath = $atom_dis_index . $_;
		s/\_.atom_atom//;
		$outpath = $out_index."/".$_.".vanatomnet";
		open(OUT,">$outpath");
		open(FILE,"$filepath");
		@atom_atoms = <FILE>;
		close(FILE);
		foreach  (@atom_atoms) {
			@atom_temp = split/\s+/,$_;
			@atom_i_ = split/\_/,$atom_temp[0];
			@atom_j_ = split/\_/,$atom_temp[1];
			$dis_ij = $atom_temp[2];
			$res_i = $atom_i_[0]; 
			$resNo_i = $atom_i_[1]; 
			$res_chain_i = $atom_i_[2]; 
			$atom_i = $atom_i_[3]; 
			$ele_i = substr($atom_i,0,1);
			$res_j = $atom_j_[0]; 
			$resNo_j = $atom_j_[1]; 
			$res_chain_j = $atom_j_[2]; 
			$atom_j = $atom_j_[3]; 
			$ele_j = substr($atom_j,0,1);
			$residue_i_site = $res_i."_".$$resNo_i;
			$residue_j_site = $res_j."_".$$resNo_j;
			$_ = $res_chain_i;
			s/a/A/;
			s/b/B/;
			s/i/I/;
			s/k/K/;
			s/n/N/;
			s/e/E/;
			s/s/S/;
			s/u/U/;
			s/z/Z/;
			$res_chain_i = $_;
			$_ = $res_chain_j;
			s/a/A/;
			s/b/B/;
			s/i/I/;
			s/k/K/;
			s/n/N/;
			s/e/E/;
			s/s/S/;
			s/u/U/;
			s/z/Z/;
			$res_chain_j = $_;
			$van_atom_dis = &van_dis($ele_i,$ele_j);
			if ( ($dis_ij <= $van_atom_dis)&&(($residue_i_site ne $residue_j_site)||($res_chain_i =~ /$res_chain_j/)) ) {
				syswrite OUT,$res_i."_".$res_chain_i."_".$resNo_i."_".$atom_i." "."van"." ";
				syswrite OUT,$res_j."_".$res_chain_j."_".$resNo_j."_".$atom_j." ".$dis_ij." ".$ele_i."\n";
			}
		}
	}
}
sub van_dis{
	my $ai = $_[0];
	my $aj = $_[1];
	my $aij = 0;
	if ( ($ai =~ /C/)&&($aj =~ /C/) ) {
		$aij = $C+$C;	}
	if ( ($ai =~ /C/)&&($aj =~ /N/) ) {
		$aij = $C+$N;   }
	if ( ($ai =~ /C/)&&($aj =~ /O/) ) {
		$aij = $C+$O;	}
	if ( ($ai =~ /C/)&&($aj =~ /S/) ) {
		$aij = $C+$S;	}
	if ( ($ai =~ /N/)&&($aj =~ /C/) ) {
		$aij = $N+$C;	}
	if ( ($ai =~ /N/)&&($aj =~ /O/) ) {
		$aij = $N+$O;	}
	if ( ($ai =~ /N/)&&($aj =~ /N/) ) {
		$aij = $N+$N;	}
	if ( ($ai =~ /N/)&&($aj =~ /S/) ) {
		$aij = $N+$S;	}
	if ( ($ai =~ /O/)&&($aj =~ /C/) ) {
		$aij = $O+$C;	}
	if ( ($ai =~ /O/)&&($aj =~ /N/) ) {
		$aij = $O+$N;	}
	if ( ($ai =~ /O/)&&($aj =~ /O/) ) {
		$aij = $O+$O;	}
	if ( ($ai =~ /O/)&&($aj =~ /S/) ) {
		$aij = $O+$S;	}
	if ( ($ai =~ /S/)&&($aj =~ /C/) ) {
		$aij = $S+$C;	}
	if ( ($ai =~ /S/)&&($aj =~ /O/) ) {
		$aij = $S+$O;	}
	if ( ($ai =~ /S/)&&($aj =~ /N/) ) {
		$aij = $S+$N;	}
	if ( ($ai =~ /S/)&&($aj =~ /S/) ) {
		$aij = $S+$S;	}
	$aij;
}
