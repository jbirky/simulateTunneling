import os
import shutil

dirs = ['wave_func/', 'wave_prob/', 'plots/phi/']

for path in dirs:
	shutil.rmtree(path)
	os.mkdir(path)