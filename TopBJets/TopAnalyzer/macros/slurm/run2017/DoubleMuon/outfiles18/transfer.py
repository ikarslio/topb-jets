import argparse, os, json, glob, time
from common import *

parser = argparse.ArgumentParser()
if __name__ == '__main__':
    parser.add_argument('--Tier2-prefix', dest='Tier2_prefix', action='store', default='root://eos.cms.rcac.purdue.edu//store/user/',help='prefix of the Tier-2 storage site')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', default=False,
             help='show verbose printouts during execution')

opts, opts_unknown = parser.parse_known_args()

log_prx = os.path.basename(__file__)+' -- '

f = open("out_files.txt", "r")
x = len(f.readlines())
f = open("out_files.txt", "r")
data = f.read()
data_list = data.split()
#print data_list
for i_tmp in range(x):
    if not os.path.isfile(data_list[i_tmp]):
        WARNING(log_prx+'temporary local output file not found: '+output_files[i_tmp])
        continue

    #TIER2_CP_CMD = 'gfal-copy -t 99 -T 99'
    TIER2_CP_CMD = 'xrdcp'
    TIER2_PREFIX = 'root://eos.cms.rcac.purdue.edu//store/user/rchawla/minitrees/run2018/DoubleMuon/'
    EXE(TIER2_CP_CMD+' '+data_list[i_tmp]+' '+TIER2_PREFIX+data_list[i_tmp])
    EXE('rm'+' '+data_list[i_tmp])
