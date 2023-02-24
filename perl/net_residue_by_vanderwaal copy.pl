#!/usr/bin/perl
my $atom_vannet_index = shift @ARGV;
my $out_index = $atom_vannet_index;
my $pdb_name = shift @ARGV;
$_ = $pdb_name;
s/\.pdb//;
@files = $_.".vanatomnet";
foreach  (@files) {
	if ( !($_ =~ /^\./) ) {
		$filepath = $atom_vannet_index."/".$_;
		s/\.vanatomnet//;
		$outpath = $out_index."/".$_.".vanresnet";
		open(OUT,">$outpath");
		open(FILE,"$filepath");
		@atoms = <FILE>;
		close(FILE);
		foreach  (@atoms) {
			@atom_temp = split/\s+/,$_;
			@atom_i_ = split/\_/,$atom_temp[0];
			@atom_j_ = split/\_/,$atom_temp[2];
			$dis_ij = $atom_temp[3];
			$res_i = $atom_i_[0]; $res_chain_i = $atom_i_[1]; $resNo_i = $atom_i_[2]; $atom_i = $atom_i_[3];
			$res_j = $atom_j_[0]; $res_chain_j = $atom_j_[1]; $resNo_j = $atom_j_[2]; $atom_j = $atom_j_[3];
			if ( abs($resNo_i - $resNo_j) == 1 ) {
				syswrite OUT, $res_i."_".$res_chain_i."_".$resNo_i." ".$atom_i."_".$resNo_i."_".$resNo_j."_".$atom_j." ".$res_j."_".$res_chain_j."_".$resNo_j." "."AA_band"." ".$dis_ij."\n";
			}
			if ( abs($resNo_i - $resNo_j) > 1 ) {
				syswrite OUT, $res_i."_".$res_chain_i."_".$resNo_i." ".$atom_i."_".$resNo_i."_".$resNo_j."_".$atom_j." ".$res_j."_".$res_chain_j."_".$resNo_j." "."AA_van"." ".$dis_ij."\n";
			}
		}
	}
}