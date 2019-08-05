import sys
sys.path.append("..")

run     = sys.argv[1]

from ac_derivation_OHC import DeriveOHC as DO

if __name__=="__main__":
    DO().generate_OHC_files(run)
    