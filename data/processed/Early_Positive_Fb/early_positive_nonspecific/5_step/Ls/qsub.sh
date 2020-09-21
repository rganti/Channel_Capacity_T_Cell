#PBS -m ae
#PBS -q short
#PBS -V
#PBS -l walltime=0:10:00,nodes=1:ppn=2 -N ODE_steps_5

cd $PBS_O_WORKDIR
echo $PBS_JOBID > job_id

python ~/SSC_python_modules/pysb_t_cell_network.py --steps 5 --early_pos_fb
