#! /usr/bin/env python
import subprocess

for i in {4}:
	subprocess.call(['make','clean'])
	subprocess.call(['make'])
	for j in {101}:
		subprocess.call(['./flow_final', 'configs/config_%d.txt' % i, 'res/conf%d_output1.txt' % i, 'res/conf%d_output2.txt' % i, 'data/data_%d.txt' % j, 'configs/model_config_%d.txt' % i]) 

