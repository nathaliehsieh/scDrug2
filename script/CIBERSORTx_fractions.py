import argparse, os, sys

## Parse command-line arguments
# process arguments
parser = argparse.ArgumentParser(description="impute cell fractions")

parser.add_argument("-i", "--input", required=True, help="path to input directory")
parser.add_argument("-o", "--output", default='./', help="path to output directory, default='./'")
parser.add_argument("-u", "--username", required=True, help="email address registered on CIBERSORTx website")
parser.add_argument("-t", "--token", required=True, help="token obtained from CIBERSORTx website")
parser.add_argument("-r", "--refsample", required=True, help="reference sample file from single cell RNA sequencing data")
parser.add_argument("-m", "--mixture", required=True, help="mixture matrix required for running CIBERSORTx")

args = parser.parse_args()


## Run CIBERSORTx fractions
os.system("docker run --rm --name cibersortx/fractions \
          -v {input_dir}:/src/data -v {output_dir}:/src/outdir \
          cibersortx/fractions --username {username} --token {token} \
          --single_cell TRUE --refsample {refsample} --mixture {mixture}".format(
          input_dir=os.path.abspath(args.input), output_dir=os.path.abspath(args.output), 
          username=args.username, token=args.token, 
          refsample=args.refsample, mixture=args.mixture))