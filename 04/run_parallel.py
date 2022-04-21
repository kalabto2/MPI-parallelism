import subprocess

FP_P_SCRIPT = "./parallel_job.sh"
ITH_LINE = 66
FP_INSTANCE_FOLDER = "./graf_mbp/"
INSTANCES_FILENAME = ["graf_10_3.txt", "graf_10_5.txt", "graf_10_6.txt", "graf_10_7.txt", "graf_12_3.txt", "graf_12_5.txt", "graf_12_6.txt", "graf_15_4.txt", "graf_15_6.txt", "graf_15_8.txt", "graf_12_9.txt"]
FP_OUT = "./out.txt"
FP_ERR = "./err.txt"

def rewrite_file (new_argument):
	f = open(FP_P_SCRIPT, "r+")
	lines = f.read().split("\n")
	lines[ITH_LINE] = 'MY_PARALLEL_PROGRAM="./a.out ' + new_argument + '"'
	f.write("\n".join(lines))
	f.close()


if __name__ == "__main__":
	# clear output_files
	with open(FP_OUT, 'r+') as f:
	    f.truncate(0)
	with open(FP_ERR, 'r+') as f:
	    f.truncate(0)

	# compile program
	subprocess.run(['mpiCC', '-fopenmp', '-std=c++20', 'main.cpp'],check=True, stdout=subprocess.PIPE, universal_newlines=True)

	# run all instances
	for instance_name in INSTANCES_FILENAME:
		rewrite_file(FP_INSTANCE_FOLDER + instance_name)
		subprocess.run(['qrun2','20c','4', 'pdp_long', FP_P_SCRIPT],check=True, stdout=subprocess.PIPE, universal_newlines=True)

