#!/usr/bin/env python

import run, multiprocessing
from copy import deepcopy as dc

states = ["al", "ak", "az", "ar", "ca", "co", "ct", "de", "dc", 
          "fl", "ga", "hi", "id", "il", "in", "ia", "ks", "ky",
          "la", "me", "md", "ma", "mi", "mn", "ms", "mo", "mt", 
	  "ne", "nv", "nh", "nj", "nm", "ny", "nc", "nd", "oh", 
	  "ok", "or", "pa", "ri", "sc", "sd", "tn", "tx", "ut", 
	  "vt", "va", "wa", "wv", "wi", "wy"]

states = ["pa", "mn", "tn", "wi", "nc", "il", "md", "tx", "wi", "va", "la", "fl"]

## POWER
args = {'no_plot': True, 'verbose': 0, 'grasp': False, 'method': 'power', 'power_restart': True,
        'borders': '', 'ncycles': 100, 'allow_trades': False, 'init': 'power:100000', 'destrand_max': 0, 
        'niter': 100, 'print_init': True, "output" : "test", 'shading': ['target'], 'nloops': 0, 'point': None, 'conv_iter': 500, 
        'tabu_length': 0, 'tol': 0.01, 'circle': '', 'ring': False, 'split_restart': False,
        'destrand_min': 2, 'destrand_inputs': False,
        "seats" : 0, "bgroup" : False, "ctol" : 0}
method = "power"

## EXCHANGE
# args = {'method': 'exchange', 'nloops': 1, 'niter': 10000, 'ncycles': 3, 'conv_iter': 500, 'ctol': 0.02, 'tol': 0.01, 
#         'shading': ['district'], 'borders': '', 'circle': '', 'ring': False, 'allow_trades': False, 'print_init': False, "output" : "test", 'no_plot': False, 
#         'bgroup': False, 'destrand_min': 5, 'point': None, 'tabu_length': 0, 'destrand_max': 50,
#         'write': '', 'seats': 0, 'power_restart': False, 'split_restart': False,
#         'grasp': False, 'destrand_inputs': False, 'verbose': 0, 'init': 'seed'}
# method = "exchange"

def run_state_dict(d): run.main(**d)

def queue_states():

  arg_list = []
  for x in range(240, 250):
    for s in states:
      arg_list.append(dc(args))
      arg_list[-1]["state"] = s
      arg_list[-1]["seed"]  = x
      arg_list[-1]["write"] = "{}/{}/s{:03d}".format(s, method, x)

  p = multiprocessing.Pool(4)
  p.map(run_state_dict, arg_list)


if __name__ == '__main__': queue_states()


