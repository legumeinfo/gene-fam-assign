#!/usr/bin/env perl

# In multicolumn input in which the first column is a mRNA ID with splice variant,
# return the data with the splice variant string stripped from this ID.
# Intended use: in tabular "gene family assignment" data derived from hmmsearch,
# produce a gene ID in the first column.
# The regex below should handle geneID.1, geneID.m1, geneID.mRNA1, geneID-T1

while (<>){
  my @F = split(/\t/, $_);
  if ( $F[0] =~ /(.+)\.\w\d+$/ ||
       $F[0] =~ /(.+\w\w+)\.\d{1,2}\.\d{1,2}$/ ||
       $F[0] =~ /(.+)\.\d+$/ ||
       $F[0] =~ /(.+)-T\d+$/ ||
       $F[0] =~ /(.+)\.mRNA\d+$/ ||
       $F[0] =~ /^([^.]+\.[^.]+\.[^.]+\.[^.]+\.[^.]+)$/ ){
         $F[0]=$1; print join("\t",@F)
       }
  else { $F[0]="XX"; print join ("\t",@F) } # print a flag; this data will need special handling
}

