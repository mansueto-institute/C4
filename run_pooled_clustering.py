#!/usr/bin/env python

import run, multiprocessing
from copy import deepcopy as dc

states = ["al", "az", "ar", "ca", "co", "ct", 
          "fl", "ga", "hi", "id", "il", "in", "ia", "ks", "ky",
          "la", "me", "md", "ma", "mi", "mn", "ms", "mo", 
	  "ne", "nv", "nh", "nj", "nm", "ny", "nc", "oh", 
	  "ok", "or", "pa", "ri", "sc", "tn", "tx", "ut", 
	  "vt", "va", "wa", "wv", "wi"]

states = ["nc"]

## POWER
power = {'no_plot': True, 'verbose': 0, 'grasp': False, 'method': 'power', 'power_restart': True,
         'borders': '', 'ncycles': 100, 'allow_trades': False, 'init': 'power:100000', 'destrand_max': 0, 
         'niter': 100, 'print_init': True, "output" : "res", 'shading': ['target'], 'nloops': 0, 'point': None, 'conv_iter': 500, 
         'tabu_length': 0, 'tol': 0.025, 'circle': '', 'ring': False, 'split_restart': False,
         'destrand_min': 2, 'destrand_inputs': False,
         "seats" : 0, "bgroup" : True, "ctol" : 0}
method = "power"

split = {'bgroup': True, 'init': 'split', 'state': 'md', 'seed': 0, 'method': 'dist_a', 'nloops': 0, 'destrand_inputs': False, 'niter': 100, 'destrand_min': 2, 'write': 'md/split/s001', 'circle': '', 'ctol': 0.02, 'power_restart': False, 'tabu_length': 0, 'ncycles': 1, 'grasp': False, 'allow_trades': False, 'borders': '', 'print_init': True, 'shading': [], 'seats': 0, 'tol': 0.025, 'no_plot': False, 'point': None, 'conv_iter': 0, 'output': 'res', 'split_restart': False, 'ring': False, 'destrand_max': 0, 'verbose': 0}


def run_state_dict(d): run.main(**d)

power["seats"] = 120
split["seats"] = 120

def queue_states():

  arg_list = []
  for x in range(271, 275):
    for s in states:
      arg_list.append(dc(power))
      arg_list[-1]["state"] = s
      arg_list[-1]["seed"]  = x
      arg_list[-1]["write"] = "{}/{}/s{:03d}".format(s + "_rep", method, x)

  p = multiprocessing.Pool(4)
  p.map(run_state_dict, arg_list)


  arg_list = []
  for s in states:
    arg_list.append(dc(split))
    arg_list[-1]["state"] = s
    arg_list[-1]["seed"]  = 1
    arg_list[-1]["write"] = "{}/split/s001".format(s + "_rep", method)

  p = multiprocessing.Pool(4)
  p.map(run_state_dict, arg_list)



if __name__ == '__main__': queue_states()


