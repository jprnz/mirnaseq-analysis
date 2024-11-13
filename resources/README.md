# Datasets
Documented here are nuances of each dataset, any modifications made, and checksums of downloaded files.

## Catalogue
Each external dataset will be descibed in `datasets.yaml` as follows

### Species codes
All miRNAs and piRNAs are expected to have a three letter species code. 
Please see the following for available species codes:
https://www.mirbase.org/help

### miRNAs 
Entries are expected to have values for `mature` and `hairpin` and each should an associated regular expression that is capable of selecting the appropriate sequences via their FASTA header. Be sure to include the `{species}` to pass along species values.
Additionally, a value for `sequence_type` is required and should be one of `dna` or `rna`.

### piRNAs
It is assumed piRNA datasets are concatenated and the same species code as can be used as the miRNA datasets.  


# Notes
## miRBase
#### Source
https://www.mirbase.org/download/
Release 22.1 - Downloaded on 08/30/2023

#### Changes
None

#### Checksums
```
922bcd22b0025877f7db2e4ed5c3571e  hairpin.fa.gz
d84a0b564217f4cdde6c54f35c017bf7  mature.fa.gz
```

## mirGeneDB
#### Source
https://www.mirgenedb.org/fasta/ALL?all=1

#### Changes
None

#### Checksums
```
6a10d9401a1cab80b32997f534b474f9  ALL.fas.gz
```

## piRNABank
#### Source
http://pirnabank.ibab.ac.in/downloads/all/
Downloaded on 08/30/23 from:

#### Changes
Changed extension to ".fa", gzipped, and concatenated

#### Checksums
```
2334c1b5385772b8d1830c6d96171eff  drosophila_pir.fa
dae94a74c4861792c70beaf444238856  human_pir.fa
b8b9bef1d5827325d47a1f431005c90b  mouse_pir.fa
95baa9d9712ea391c644070f9cb4bc70  rat_pir.fa
```

## piRNAdb
#### Source
https://www.pirnadb.org/download/archive
Downloaded on 08/30/2023

#### Changes
File were concatenated

#### Checksums
```
d137de21b2bedde67f324e8b9aa4cd73  piRNAdb.cel.v1_7_6.fa.gz
834609d656f9b029ec47df9982d4d798  piRNAdb.cgr.v1_7_6.fa.gz
8574f7e46c515bbf178d4bd15a5d47c2  piRNAdb.dme.v1_7_6.fa.gz
5b706d61aa21c5ff83818375d1070d1b  piRNAdb.hsa.v1_7_6.fa.gz
436a2e27955fb167b2dd4a5dc768285b  piRNAdb.mmu.v1_7_6.fa.gz
825326bb8e04e7d8745615b51d8b6dc6  piRNAdb.rno.v1_7_6.fa.gz
```
