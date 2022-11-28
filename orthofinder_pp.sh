#PBS -N Orthofinder
#PBS -l nodes=1:ppn=8,vmem=64gb,walltime=800:00:00
#PBS -q default
#PBS -V
#PBS -o Ortho.out
#PBS -e Ortho.err

# modulos
module load orthofinder/2.4.0
module load mcl/14-137
module load ncbi-blast+/2.9.0
module load mafft/7.305b

# ubicacion
cd $PBS_O_WORKDIR

orthofinder -f primary_transcripts/ -t 8 -a 8 -S blast -I 1.7 -M msa -A mafft -T fasttree

