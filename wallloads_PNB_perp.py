import os
import plot_walloads 

dir='/home/vallar/'
if os.uname().nodename!='spcpc182':
    dir='/home/matval/'
dir+='WORK/ASCOT/runs/SA_003/ripple/perp/runs_072020/'
fn=f'{dir}/ascot_towall_output_allwall.h5'
flag_fild=11

walloads_on_fild,x1x2x3, y1y2y3, z1z2z3 = plot_wallloads.loads_on_flag(fn,flag_fild)