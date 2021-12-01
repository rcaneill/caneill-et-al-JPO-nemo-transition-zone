from pathlib import Path
from os import chdir, system
from shutil import copy, copytree, move
import csv
import f90nml

def good_type(s):
    """
    Converts the string s to either a float, a bool or a string
    """
    if '.true.' in s:
        return True
    elif '.false.' in s:
        return False
    try:
        return int(s)
    except ValueError:
        try:
            return float(s)
        except ValueError:
            return s
        
            

if __name__ == '__main__':
    with open(snakemake.input['cfg'], 'r') as file:
        reader = csv.reader(file)
        headers = next(reader)
        headers = [i.strip() for i in headers]
        for _cfg in reader:
            cfg = [i.strip() for i in _cfg]
            # Create the config name
            name = 'EXP'
            for (nam, value) in zip(headers[1:], cfg[1:]):
                name += f'_{nam}_{value}'
                
            # Create the config folder from EXP00
            exp_dir = Path(snakemake.params['exp']) / 'Experiments'
            copytree(Path(snakemake.params['exp']) / 'EXP00', exp_dir / name)
                        
            # Add config to experiments.csv
            with open(snakemake.input['expcsv'], 'a') as f:
                f.write(f'{name}, {cfg[0]}\n')

            # Change namelist
            ###chdir(exp_dir / name)
            namelist = f90nml.read(exp_dir / name / 'namelist_cfg')
            # We search for the name of the sub-namelist for each param
            # Not efficient method but sufficiently good for us
            sn = {}
            for i in namelist.keys():
                sn.update(dict.fromkeys(namelist[i], i))

            patch = {}
            # We create a patch to keep the namelist formatting
            for (nam, value) in zip(headers[1:], cfg[1:]):
                if sn[nam] not in patch.keys():
                    patch[sn[nam]] = {}
                patch[sn[nam]][nam] = good_type(value)
            
            # Choose timestep, executing time, number of processors, etc
            try:
                res = patch['namusr_def']['rn_e1_deg']
            except KeyError:
                res = namelist['namusr_def']['rn_e1_deg'] # 1 deg or 2 deg
            if res == 1.:
                # namelist
                rn_rdt = 2700.
                nn_iglo = 42
                nn_jglo = 79
                rn_Ld = 100.e+3
                rn_Le = 100.e+3
                rn_Lv = 100.e+3
                # run_nemo
                TIME = '3-00:00:00'
                N_PROCS_TOT = 32
                SBATCH_PROJ = snakemake.config['sbatch_proj']
                N_TIME_STEP_PER_YEAR = 11520
                N_PROCS_NEMO = 28
                N_PROCS_XIOS = 4
            elif res == 2.:
                # namelist
                rn_rdt = 2700. * 2
                nn_iglo = 20
                nn_jglo = 40
                rn_Ld = 200.e+3
                rn_Le = 200.e+3
                rn_Lv = 200.e+3
                # run_nemo
                TIME = '30:00:00'
                N_PROCS_TOT = 9
                SBATCH_PROJ = snakemake.config['sbatch_proj']
                N_TIME_STEP_PER_YEAR = 5760
                N_PROCS_NEMO = 8
                N_PROCS_XIOS = 1
            else:
                raise(ValueError('Only 1 and 2 degrees resolution are authorized'))
            # change patch
            for i in ['namdom', 'namusr_def', 'namtra_ldf', 'namtra_eiv', 'namdyn_ldf']:
                if not i in patch.keys():
                    patch[i] = {}
            patch[sn['rn_rdt']]['rn_rdt'] = rn_rdt
            patch[sn['nn_iglo']]['nn_iglo'] = nn_iglo
            patch[sn['nn_jglo']]['nn_jglo'] = nn_jglo
            patch[sn['rn_ld']]['rn_ld'] = rn_Ld
            patch[sn['rn_le']]['rn_le'] = rn_Le
            patch[sn['rn_lv']]['rn_lv'] = rn_Lv
            f90nml.patch(
                Path(snakemake.params['exp']) / 'EXP00/namelist_cfg',
                patch,
                exp_dir / name / 'namelist_cfg'
            )
            
            # Change run_nemo_restart.sh
            with open(exp_dir / name / 'run_nemo_restarts_template.sh', 'r') as f:
                run_nemo = f.read()
            # TODO to finish
            run_nemo = run_nemo.format(
                TIME=TIME,
                N_PROCS_TOT=N_PROCS_TOT,
                SBATCH_PROJ=SBATCH_PROJ,
                N_TIME_STEP_PER_YEAR=N_TIME_STEP_PER_YEAR,
                N_PROCS_NEMO=N_PROCS_NEMO,
                N_PROCS_XIOS=N_PROCS_XIOS
            )
            with open(exp_dir / name / 'run_nemo_restart.sh', 'w') as f:
                f.write(run_nemo)
            
            
