def preparation(conf, clean_mode = 'grep'):
    name, ntraj = conf.PDB_filename, conf.ntraj
    import os
    current_dir = os.getcwd()
    if not os.path.exists(conf.jobname):
        os.mkdir(conf.jobname)
    os.chdir(conf.jobname)
    pdb_name = name[:-4]
    
    match clean_mode:
        case 'grep':
            os.system(f'''grep "^ATOM" ../{pdb_name}.pdb > {pdb_name}_clean.pdb''')
        case 'pyrosetta':
            import pyrosetta
            pyrosetta.init("-mute all")
            pyrosetta.toolbox.cleanATOM(f"../{conf.PDB_filename}")
            os.rename(f'../{pdb_name}.clean.pdb', f'{pdb_name}_clean.pdb')
        case 'custom':
            from script.clean_pdb import main as clean_pdb
            clean_pdb(f'../{pdb_name}.pdb')
            os.rename(f'../{conf.PDB_filename}.clean.pdb', f'{pdb_name}_clean.pdb')
    
    for i in range(1, ntraj + 1):
        if not os.path.exists(str(i)):
            os.mkdir(str(i))
            os.system(f'touch {str(i)}/stdout.txt')
    os.chdir(current_dir)
    return 0
def relax_job(args):
    name_clean, nstruct, scorefxn_name, dest = args
    import os
    import logging
    current_dir = os.getcwd()
    os.chdir(dest)
    logging.basicConfig(filename='stdout.txt', level=logging.INFO)
    import pyrosetta
    pyrosetta.init()
    pose = pyrosetta.pose_from_pdb('../' + name_clean)
    scorefxn = pyrosetta.create_score_function(scorefxn_name)
    xml = pyrosetta.rosetta.protocols.rosetta_scripts.XmlObjects.create_from_string("""
    <ROSETTASCRIPTS>
        <SCOREFXNS>
        <ScoreFunction name="SFX1" weights="ref2015_cart">
           <Reweight scoretype="coordinate_constraint" weight="1.0"/>
        </ScoreFunction>
        </SCOREFXNS>
        <RESIDUE_SELECTORS>
        </RESIDUE_SELECTORS>
        <TASKOPERATIONS>
        </TASKOPERATIONS>
        <FILTERS>
        </FILTERS>
        <MOVERS>
           <AtomCoordinateCstMover name="coord_cst" />
           <FastRelax name="relax" cartesian="true" scorefxn="SFX1" />
        </MOVERS>
        <APPLY_TO_POSE/>
        <PROTOCOLS>
           <Add mover="coord_cst" />
           <Add mover="relax" />
        </PROTOCOLS>
    </ROSETTASCRIPTS>
    """).get_mover("ParsedProtocol")
    
    working_dir = os.getcwd()
    output_dir = dest
    jd = pyrosetta.toolbox.py_jobdistributor.PyJobDistributor(pdb_name=name_clean[:-4], nstruct=nstruct, scorefxn=scorefxn)
    jd.native_pose = pose
    while not jd.job_complete:
        test_pose = pose.clone()
        xml.apply(test_pose)
        jd.output_decoy(test_pose)
    os.chdir(current_dir)
    return 0
def run_relax_job_parallel(conf):
    name, nstruct, scorefxn_name, num_th = conf.PDB_filename, conf.nstruct, conf.relax_scorefxn, conf.num_threads
    name = name[:-4] + '_clean.pdb'
    ntraj = conf.ntraj
    if not conf.coord_cst:
        conf.jobname = conf.jobname + '_without_cst'
        preparation(conf)
    from multiprocessing.pool import Pool
    import os
    os.chdir(conf.jobname)
    args = [(name, nstruct, scorefxn_name, str(i)) for i in range(1, ntraj + 1)]
    with Pool(num_th) as pool:
        pool.map(relax_job, args, chunksize=1)
    os.chdir('../')
    if not conf.coord_cst:
        conf.jobname = conf.jobname[:-12]
    return 0

def get_scores(conf):
    conf.back()
    import pandas as pd
    import glob
    import os
    os.chdir(conf.jobname)
    dall_scores = pd.DataFrame(columns = ['description', 'total_score'])
    for f in glob.glob('*/*.fasc'):
        pth = str(f).split('/')[0]
        with open(f, 'r') as file:
            tmp = pd.DataFrame()
            for line in file.readlines():
                total_score = float(line[line.find('total_score'):line.find('"yhh_planarity')].split(':')[1].strip()[:-1])
                name_decoy = line[line.find("decoy"):line.find(', "filename"')].split(':')[-1].split('"')[1]
                dic = {'description': [pth + '/' + name_decoy], 'total_score':[total_score]}
                tmp = pd.DataFrame(dic)
                dall_scores = pd.concat([dall_scores, tmp], ignore_index = True)

    dres = dall_scores[dall_scores.total_score == dall_scores.total_score.min()]
    source  = list(dres.description)[0]
    pdb = source.split('/')[1]
    if not os.path.exists('mutations'):
        os.mkdir('mutations')
    target = "mutations/relaxed_"  + pdb
    os.system(f'cp {source} {target}')
    os.chdir('../')
    return target.split('/')[-1], dall_scores.total_score.min()
def mutfiles_preparation(conf):
    import os
    conf.back()
    os.chdir(f'{conf.jobname}/mutations')
    mut = conf.mut
    aa_list = conf.aa_list
    mutations = []
    for aa in range(0, len(mut)):
        pos = mut[aa]
        AA = mut[aa][0]
        if not os.path.exists(pos):
            os.makedirs(pos)
        for m in aa_list:
            if m != AA:
                mutation = pos+m
                mutations.append(mutation)
                if not os.path.exists(pos+"/"+mutation):
                    os.makedirs(pos+"/"+mutation)
                with open(pos+"/"+mutation+"/mutfile",'w') as mut_file:
                    mut_file.write("total 1\n")
                    mut_file.write("1\n")
                    mut_file.write("%s %d %s\n" %(AA, int(mut[aa][1:]), m))
    os.chdir('../../')
    return 0

def ddg_job(args):
    pdb_name, py_flags, dest = args
    import os
    import logging
    current_dir = os.getcwd()
    os.chdir(dest)
    logging.basicConfig(filename='stdout.txt', level=logging.INFO)
    import pyrosetta
    
    pyrosetta.init(py_flags)
    pose = pyrosetta.pose_from_pdb('../../' + pdb_name)
    pyrosetta.rosetta.protocols.ddg.CartesianddG.run(pose)
    os.chdir(current_dir)
    return 0

def run_ddg_calc_parallel(conf):
    conf.back()
    import os 
    os.chdir(f'{conf.jobname}/mutations')
    name = conf.cleanaxed_pdb
    num_th = conf.num_threads
    import glob
    args = [(name, conf.py_flags, x[:-7]) for x in glob.glob('*/*/mutfile')]
    from multiprocessing.pool import Pool
    with Pool(num_th) as pool:
        pool.map(ddg_job, args, chunksize=1)
    os.chdir('../../')
    return 0
def analyze_ssm(conf):
    conf.back()
    mut = conf.mut
    aa_list = conf.aa_list
    import os
    import glob
    import numpy as np
    os.chdir(f'{conf.jobname}/mutations')
    def nice_order(aa_l):
        nice_order_for_heatmap  = ["G","P","E","D","R","K","H","Q","N","T","S","Y","W","F","M","C","I","L","V","A"]
        aa_list_for_an = nice_order_for_heatmap.copy()
        for item in nice_order_for_heatmap:
            if item not in aa_l:
                aa_list_for_an.remove(item)
        return aa_list_for_an
    aa_list = nice_order(aa_list)
    ala_scan = {}
    ssm = np.zeros([len(aa_list), len(mut)], dtype=float)
    pdb_name = conf.cleanaxed_pdb[:-4]

    for aa in range(0, len(mut)):
        pos = mut[aa]
        AA = mut[aa][0]
        for i in range(0, len(aa_list)):
            amino_acid = aa_list[i]
            if AA != amino_acid:
                mutation = mut[aa]+amino_acid
                n_WT = 0
                n_MUT = 0
                score_WT = 0
                score_MUT = 0
                with open(pos+"/"+mutation+"/mutfile.ddg", 'r') as mutfile:
                    for line in mutfile:
                        if "WT" in line:
                            score = float(line.split()[3])
                            n_WT += 1
                            score_WT += score
                        elif "MUT_" in line:
                            score = float(line.split()[3])
                            n_MUT += 1
                            score_MUT += score
                score_WT = score_WT/n_WT
                score_MUT = score_MUT/n_MUT
                ddG = score_MUT - score_WT
                ssm[i,aa] = float(ddG)


    np.savetxt(f"../../{conf.jobname}_SSM_ddg.csv", ssm, delimiter=",")
    os.chdir('../../')
    return 0



def plot_heatmap(conf):
    
    import numpy as np
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt 
    import os
    conf.back()
    def nice_order(aa_l):
        nice_order_for_heatmap  = ["G","P","E","D","R","K","H","Q","N","T","S","Y","W","F","M","C","I","L","V","A"]
        aa_list_for_an = nice_order_for_heatmap.copy()
        for item in nice_order_for_heatmap:
            if item not in aa_l:
                aa_list_for_an.remove(item)
        return aa_list_for_an
    aa_list = nice_order(conf.aa_list)
    mut = conf.mut
    df_ddg = pd.read_csv(f'{conf.jobname}_SSM_ddg.csv', header=None)    
    rows = {i:aa_list[i] for i in range(len(aa_list))}
    columns = {i:mut[i] for i in range(len(mut))}
    
    df_ddg = df_ddg.rename(columns=columns)
    df_ddg = df_ddg.rename(index=rows)
    minc = df_ddg.min().min()
    range_color = [minc, (-1.5) * minc]
    sns.set (rc = {'figure.figsize':(15, 10)})
    ax = sns.heatmap(df_ddg,  cmap="jet", annot=True, fmt=".1f",linewidth=.5, vmin = minc, vmax = (-1.5) * minc, \
                     cbar_kws={'label': 'ΔΔG', 'orientation': 'vertical'})#, annot_kws={"size": 20}
    ax.set(xlabel="Mutations", ylabel="Amino Acid")
    ax.xaxis.tick_top()

