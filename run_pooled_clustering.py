#!/usr/bin/env python

import run, multiprocessing
from copy import deepcopy as dc

states = ["al", "ak", "az", "ar", "ca", "co", "ct", "de", "dc", 
          "fl", "ga", "hi", "id", "il", "in", "ia", "ks", "ky",
          "la", "me", "md", "ma", "mi", "mn", "ms", "mo", "mt", 
	  "ne", "nv", "nh", "nj", "nm", "ny", "nc", "nd", "oh", 
	  "ok", "or", "pa", "ri", "sc", "sd", "tn", "tx", "ut", 
	  "vt", "va", "wa", "wv", "wi", "wy"]

states = ["ca", "az", "fl", "pa", "il"]

args = {'no_plot': False, 'verbose': 0, 'grasp': False, 'method': 'power', 'power_restart': True,
        'borders': '', 'ncycles': 100, 'allow_trades': False, 'init': 'power:100000', 'destrand_max': 0, 
        'niter': 100, 'print_init': True, 'shading': ['target'], 'nloops': 0, 'point': None, 'conv_iter': 500, 
        'tabu_length': 0, 'tol': 0.01, 'circle': '', 'ring': False, 'split_restart': False,
        'destrand_min': 2, 'destrand_inputs': False,
        "seats" : 0, "blocks" : False, "ctol" : 0}

def run_state_dict(d): run.main(**d)

def queue_states():

  arg_list = []
  for x in [220]:
    for s in states:
      arg_list.append(dc(args))
      arg_list[-1]["state"] = s
      arg_list[-1]["seed"]  = x
      arg_list[-1]["write"] = "{}/power/s{:03d}".format(s, x)

  p = multiprocessing.Pool(4)
  p.map(run_state_dict, arg_list)


if __name__ == '__main__': queue_states()


