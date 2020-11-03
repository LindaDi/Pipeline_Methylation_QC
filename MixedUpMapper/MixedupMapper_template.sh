# Linda Dieckmann, 10/2020
# This script is an attempt to make it easier for you to apply mixup mapper when following the QC methylation pipeline.
# Please always cite the original software used.

# for further information take a look at: 
https://github.com/molgenis/systemsgenetics/wiki/Resolving-mixups


###############################
#######prepare genotypes#######
##############################

#use already prepared genotype file
cp file/path/to/genotypes.* file/path/to/MixUpfolder

#use GenotypeHarmonizer to convert to TriTyper format
java -jar MixedUpMapper/GenotypeHarmonizer-1.4.20-SNAPSHOT/GenotypeHarmonizer.jar -i file/path/to/genotypes -I PLINK_BED -o /file/path/to/folder/harmonizer/in/your/folder -O TRITYPER  --update-id


##########################################
#########prepare methylation data#########
##########################################
#use own normalized/batch-corrected data
#center CpG-probes
#attention: format for methylation data: IDs in columns, probes in rows, first column-name is an empty tab!

java -jar MixedUpMapper/eqtl-mapping-pipeline-1.2.4E-SNAPSHOT/eqtl-mapping-pipeline.jar --mode normalize --in methylation/txt/prepared/for/Mixup  --out finalbeta_norm --centerscale

##annotation file is tab-sep file with header and CpGs in rows (Platform, Ill450k_ArrayAddress, Chr, ChrStart,ChrEnd)
##genotypemethylationcoupling: genotypeIDs and methylation IDs (two columns, no header)

##########################
####Mixedupmapper#######
#########################
java -Xmx15g -Xms15g -jar MixedUpMapper/eqtl-mapping-pipeline-1.2.4E-SNAPSHOT/eqtl-mapping-pipeline.jar --mode mixupmapper --in folder/harmonizer --out mixup --inexp txt.gz/file/you/gotout/in/finalbeta_norm/before --inexpplatform EPIC --inexpannot annotation.txt --gte genotypemethylationcoupling.txt


###check file mixup/BestMatchPerGenotype.txt###
grep true BestMatchPerGenotype.txt

## according to genotype, which IDs in methylation data fit best (from their methylation profile, meQTLs) to the respective genotype

grep true BestMatchPerTrait.txt





