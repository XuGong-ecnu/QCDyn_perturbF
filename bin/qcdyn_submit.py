#!/usr/bin/env python

"""
Created on Sun. Mar. 12 2021 @NYU-SH
Updated on Tue. Nov. 05 2021 @NYU-SH

@author: Zhubin Hu
"""
import subprocess
import argparse
import shutil
import os

parser = argparse.ArgumentParser(description='Submit the single or multi QCDyn jobs to NYU-SH HPC')
parser.add_argument('-p', '--part',help='specify the partition to submit, argon/parallel/serial, default is argon',type=str,default='argon')
parser.add_argument('-w', '--node',help='specify the particular node to submit, default is None',type=str,default=None)
parser.add_argument('-c', '--constraint',help='specify the specific cpu model to submit, default is None',type=str,default=None)
parser.add_argument('-n','--nprocs',help='specify the number of jobs, when n > 1, multi jobs will be submited, default is 1 (single job)',nargs='?',default='1',const='1')
parser.add_argument('-m','--mem',help='specify the size of memory for each job, default is 5, unit is GB',nargs='?',default='5',const='5')
parser.add_argument('-g','--gpu',help='set whether to submit a GPU job or not, default is False',default=False,action='store_true')
parser.add_argument('-o','--options',help='specify options of QCDyn, which will override the parameters in input control file, e.g., "--nsteps=100", default is None',type=str,default=None)
parser.add_argument('input_file',nargs='*',help='input the name of QCDyn input control files (*.inp)')
args = parser.parse_args()

def sh_script_gen(part,node,constraint,nprocs,mem,gpu,options,input_file,count):
    base_name = os.path.splitext(input_file)[0]
    if (nprocs != '1'):
        base_name = base_name + '_' + str(count)
    script_name = base_name + '.sh'
    with open(script_name,'w') as sh:
        sh.write('#!/bin/bash\n')
        sh.write('#SBATCH --job-name=' + 'QCDyn-' + base_name + '\n')
        sh.write('#SBATCH --output=' + base_name + '.%j.o' + '\n')
        sh.write('#SBATCH --error=' + base_name + '.%j.e' + '\n')
        sh.write('#SBATCH -p ' + part + '\n')
        if node != None:
            sh.write('#SBATCH -w ' + node + '\n')
        if constraint != None:
            sh.write('#SBATCH --constraint=' + constraint + '\n')
        if part == 'parallel':
            sh.write('#SBATCH --time=20-00:00:00\n')
        if part == 'serial':
            sh.write('#SBATCH --time=10-00:00:00\n')
        if part == 'debug':
            sh.write('#SBATCH --time=1-00:00:00\n')
        if part == 'argon':
            sh.write('#SBATCH --qos=argon\n')
        if gpu: # used for GPU job only
            sh.write('#SBATCH --gres=gpu:1\n')
        sh.write('#SBATCH --nodes=1\n')
        sh.write('#SBATCH -c 1\n')
        sh.write('#SBATCH --mem=' + mem + 'GB' + '\n\n')
        sh.write('export MODULEPATH=$MODULEPATH:/xspace/sungroup/modules/\n')
        sh.write('module purge\n')
        sh.write('module use /xspace/zl3289/modules/\n')
        sh.write('module load qcdyn/0.3.dev.cesare\n')
        if options != None:
            sh.write('QCDyn ' + options + ' ' + input_file + ' > ' + base_name + '.log' + '\n')
        else:
            sh.write('QCDyn ' + input_file + ' > ' + base_name + '.log' + '\n')
    return script_name

for i in range(len(args.input_file)):
    if (args.nprocs == '1'): # Single jobs
        script_name = sh_script_gen(args.part,args.node,args.constraint,args.nprocs,args.mem,args.gpu,args.options,args.input_file[i],0)
        sub = "sbatch  " + script_name
        subp = subprocess.Popen(sub,shell=True)
        subp.wait()
    else: # Multi-jobs using same input control file
        base_name = args.input_file[i][:-4]
        curdir = os.getcwd()
        tmpdir = base_name + '_jobs'
        if (os.path.exists(tmpdir)): # if old dir exists, rename it
            os.rename(tmpdir, tmpdir+'_bak')
        os.mkdir(tmpdir)
        os.chdir(tmpdir)
        for n in range(1, int(args.nprocs)+1):
            os.mkdir(str(n))
            os.chdir(str(n))
            shutil.copyfile(curdir + '/' + args.input_file[i], './' + args.input_file[i])
            script_name = sh_script_gen(args.part,args.node,args.constraint,args.nprocs,args.mem,args.gpu,args.options,args.input_file[i],n)
            sub = "sbatch  " + script_name
            subp = subprocess.Popen(sub,shell=True)
            subp.wait()
            os.chdir('../') # back to _jobs dir
        os.chdir('../')     # back to directory of input file
